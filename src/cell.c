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
#include "active.h"
#include "atomic.h"
#include "chemistry.h"
#include "drift.h"
#include "engine.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "memswap.h"
#include "minmax.h"
#include "scheduler.h"
#include "space.h"
#include "space_getsid.h"
#include "stars.h"
#include "timers.h"
#include "tools.h"
#include "tracers.h"

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
 * @brief Link the cells recursively to the given #part array.
 *
 * @param c The #cell.
 * @param parts The #part array.
 *
 * @return The number of particles linked.
 */
int cell_link_parts(struct cell *c, struct part *parts) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank)
    error("Linking foreign particles in a local cell!");

  if (c->hydro.parts != NULL)
    error("Linking parts into a cell that was already linked");
#endif

  c->hydro.parts = parts;

  /* Fill the progeny recursively, depth-first. */
  if (c->split) {
    int offset = 0;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        offset += cell_link_parts(c->progeny[k], &parts[offset]);
    }
  }

  /* Return the total number of linked particles. */
  return c->hydro.count;
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

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank)
    error("Linking foreign particles in a local cell!");

  if (c->grav.parts != NULL)
    error("Linking gparts into a cell that was already linked");
#endif

  c->grav.parts = gparts;

  /* Fill the progeny recursively, depth-first. */
  if (c->split) {
    int offset = 0;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        offset += cell_link_gparts(c->progeny[k], &gparts[offset]);
    }
  }

  /* Return the total number of linked particles. */
  return c->grav.count;
}

/**
 * @brief Link the cells recursively to the given #spart array.
 *
 * @param c The #cell.
 * @param sparts The #spart array.
 *
 * @return The number of particles linked.
 */
int cell_link_sparts(struct cell *c, struct spart *sparts) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank)
    error("Linking foreign particles in a local cell!");

  if (c->stars.parts != NULL)
    error("Linking sparts into a cell that was already linked");
#endif

  c->stars.parts = sparts;

  /* Fill the progeny recursively, depth-first. */
  if (c->split) {
    int offset = 0;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        offset += cell_link_sparts(c->progeny[k], &sparts[offset]);
    }
  }

  /* Return the total number of linked particles. */
  return c->stars.count;
}

/**
 * @brief Recurse down foreign cells until reaching one with hydro
 * tasks; then trigger the linking of the #part array from that
 * level.
 *
 * @param c The #cell.
 * @param parts The #part array.
 *
 * @return The number of particles linked.
 */
int cell_link_foreign_parts(struct cell *c, struct part *parts) {

#ifdef WITH_MPI

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank)
    error("Linking foreign particles in a local cell!");
#endif

  /* Do we have a hydro task at this level? */
  if (c->mpi.hydro.recv_xv != NULL) {

    /* Recursively attach the parts */
    const int counts = cell_link_parts(c, parts);
#ifdef SWIFT_DEBUG_CHECKS
    if (counts != c->hydro.count)
      error("Something is wrong with the foreign counts");
#endif
    return counts;
  }

  /* Go deeper to find the level where the tasks are */
  if (c->split) {
    int count = 0;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        count += cell_link_foreign_parts(c->progeny[k], &parts[count]);
      }
    }
    return count;
  } else {
    return 0;
  }

#else
  error("Calling linking of foregin particles in non-MPI mode.");
#endif
}

/**
 * @brief Recurse down foreign cells until reaching one with gravity
 * tasks; then trigger the linking of the #gpart array from that
 * level.
 *
 * @param c The #cell.
 * @param gparts The #gpart array.
 *
 * @return The number of particles linked.
 */
int cell_link_foreign_gparts(struct cell *c, struct gpart *gparts) {

#ifdef WITH_MPI

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank)
    error("Linking foreign particles in a local cell!");
#endif

  /* Do we have a hydro task at this level? */
  if (c->mpi.grav.recv != NULL) {

    /* Recursively attach the gparts */
    const int counts = cell_link_gparts(c, gparts);
#ifdef SWIFT_DEBUG_CHECKS
    if (counts != c->grav.count)
      error("Something is wrong with the foreign counts");
#endif
    return counts;
  }

  /* Go deeper to find the level where the tasks are */
  if (c->split) {
    int count = 0;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        count += cell_link_foreign_gparts(c->progeny[k], &gparts[count]);
      }
    }
    return count;
  } else {
    return 0;
  }

#else
  error("Calling linking of foregin particles in non-MPI mode.");
#endif
}

/**
 * @brief Recursively count the number of #part in foreign cells that
 * are in cells with hydro-related tasks.
 *
 * @param c The #cell.
 *
 * @return The number of particles linked.
 */
int cell_count_parts_for_tasks(const struct cell *c) {

#ifdef WITH_MPI

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank)
    error("Counting foreign particles in a local cell!");
#endif

  /* Do we have a hydro task at this level? */
  if (c->mpi.hydro.recv_xv != NULL) {
    return c->hydro.count;
  }

  if (c->split) {
    int count = 0;
    for (int k = 0; k < 8; ++k) {
      if (c->progeny[k] != NULL) {
        count += cell_count_parts_for_tasks(c->progeny[k]);
      }
    }
    return count;
  } else {
    return 0;
  }

#else
  error("Calling linking of foregin particles in non-MPI mode.");
#endif
}

/**
 * @brief Recursively count the number of #gpart in foreign cells that
 * are in cells with gravity-related tasks.
 *
 * @param c The #cell.
 *
 * @return The number of particles linked.
 */
int cell_count_gparts_for_tasks(const struct cell *c) {

#ifdef WITH_MPI

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank)
    error("Counting foreign particles in a local cell!");
#endif

  /* Do we have a hydro task at this level? */
  if (c->mpi.grav.recv != NULL) {
    return c->grav.count;
  }

  if (c->split) {
    int count = 0;
    for (int k = 0; k < 8; ++k) {
      if (c->progeny[k] != NULL) {
        count += cell_count_gparts_for_tasks(c->progeny[k]);
      }
    }
    return count;
  } else {
    return 0;
  }

#else
  error("Calling linking of foregin particles in non-MPI mode.");
#endif
}

/**
 * @brief Pack the data of the given cell and all it's sub-cells.
 *
 * @param c The #cell.
 * @param pc Pointer to an array of packed cells in which the
 *      cells will be packed.
 * @param with_gravity Are we running with gravity and hence need
 *      to exchange multipoles?
 *
 * @return The number of packed cells.
 */
int cell_pack(struct cell *restrict c, struct pcell *restrict pc,
              const int with_gravity) {

#ifdef WITH_MPI

  /* Start by packing the data of the current cell. */
  pc->hydro.h_max = c->hydro.h_max;
  pc->stars.h_max = c->stars.h_max;
  pc->hydro.ti_end_min = c->hydro.ti_end_min;
  pc->hydro.ti_end_max = c->hydro.ti_end_max;
  pc->grav.ti_end_min = c->grav.ti_end_min;
  pc->grav.ti_end_max = c->grav.ti_end_max;
  pc->stars.ti_end_min = c->stars.ti_end_min;
  pc->stars.ti_end_max = c->stars.ti_end_max;
  pc->hydro.ti_old_part = c->hydro.ti_old_part;
  pc->grav.ti_old_part = c->grav.ti_old_part;
  pc->grav.ti_old_multipole = c->grav.ti_old_multipole;
  pc->stars.ti_old_part = c->stars.ti_old_part;
  pc->hydro.count = c->hydro.count;
  pc->grav.count = c->grav.count;
  pc->stars.count = c->stars.count;
  pc->maxdepth = c->maxdepth;

  /* Copy the Multipole related information */
  if (with_gravity) {
    const struct gravity_tensors *mp = c->grav.multipole;

    pc->grav.m_pole = mp->m_pole;
    pc->grav.CoM[0] = mp->CoM[0];
    pc->grav.CoM[1] = mp->CoM[1];
    pc->grav.CoM[2] = mp->CoM[2];
    pc->grav.CoM_rebuild[0] = mp->CoM_rebuild[0];
    pc->grav.CoM_rebuild[1] = mp->CoM_rebuild[1];
    pc->grav.CoM_rebuild[2] = mp->CoM_rebuild[2];
    pc->grav.r_max = mp->r_max;
    pc->grav.r_max_rebuild = mp->r_max_rebuild;
  }

#ifdef SWIFT_DEBUG_CHECKS
  pc->cellID = c->cellID;
#endif

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      pc->progeny[k] = count;
      count += cell_pack(c->progeny[k], &pc[count], with_gravity);
    } else {
      pc->progeny[k] = -1;
    }

  /* Return the number of packed cells used. */
  c->mpi.pcell_size = count;
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Pack the tag of the given cell and all it's sub-cells.
 *
 * @param c The #cell.
 * @param tags Pointer to an array of packed tags.
 *
 * @return The number of packed tags.
 */
int cell_pack_tags(const struct cell *c, int *tags) {

#ifdef WITH_MPI

  /* Start by packing the data of the current cell. */
  tags[0] = c->mpi.tag;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL)
      count += cell_pack_tags(c->progeny[k], &tags[count]);

#ifdef SWIFT_DEBUG_CHECKS
  if (c->mpi.pcell_size != count) error("Inconsistent tag and pcell count!");
#endif  // SWIFT_DEBUG_CHECKS

  /* Return the number of packed tags used. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Unpack the data of a given cell and its sub-cells.
 *
 * @param pc An array of packed #pcell.
 * @param c The #cell in which to unpack the #pcell.
 * @param s The #space in which the cells are created.
 * @param with_gravity Are we running with gravity and hence need
 *      to exchange multipoles?
 *
 * @return The number of cells created.
 */
int cell_unpack(struct pcell *restrict pc, struct cell *restrict c,
                struct space *restrict s, const int with_gravity) {

#ifdef WITH_MPI

  /* Unpack the current pcell. */
  c->hydro.h_max = pc->hydro.h_max;
  c->stars.h_max = pc->stars.h_max;
  c->hydro.ti_end_min = pc->hydro.ti_end_min;
  c->hydro.ti_end_max = pc->hydro.ti_end_max;
  c->grav.ti_end_min = pc->grav.ti_end_min;
  c->grav.ti_end_max = pc->grav.ti_end_max;
  c->stars.ti_end_min = pc->stars.ti_end_min;
  c->stars.ti_end_max = pc->stars.ti_end_max;
  c->hydro.ti_old_part = pc->hydro.ti_old_part;
  c->grav.ti_old_part = pc->grav.ti_old_part;
  c->grav.ti_old_multipole = pc->grav.ti_old_multipole;
  c->stars.ti_old_part = pc->stars.ti_old_part;
  c->hydro.count = pc->hydro.count;
  c->grav.count = pc->grav.count;
  c->stars.count = pc->stars.count;
  c->maxdepth = pc->maxdepth;

#ifdef SWIFT_DEBUG_CHECKS
  c->cellID = pc->cellID;
#endif

  /* Copy the Multipole related information */
  if (with_gravity) {

    struct gravity_tensors *mp = c->grav.multipole;

    mp->m_pole = pc->grav.m_pole;
    mp->CoM[0] = pc->grav.CoM[0];
    mp->CoM[1] = pc->grav.CoM[1];
    mp->CoM[2] = pc->grav.CoM[2];
    mp->CoM_rebuild[0] = pc->grav.CoM_rebuild[0];
    mp->CoM_rebuild[1] = pc->grav.CoM_rebuild[1];
    mp->CoM_rebuild[2] = pc->grav.CoM_rebuild[2];
    mp->r_max = pc->grav.r_max;
    mp->r_max_rebuild = pc->grav.r_max_rebuild;
  }

  /* Number of new cells created. */
  int count = 1;

  /* Fill the progeny recursively, depth-first. */
  c->split = 0;
  for (int k = 0; k < 8; k++)
    if (pc->progeny[k] >= 0) {
      struct cell *temp;
      space_getcells(s, 1, &temp);
      temp->hydro.count = 0;
      temp->grav.count = 0;
      temp->stars.count = 0;
      temp->loc[0] = c->loc[0];
      temp->loc[1] = c->loc[1];
      temp->loc[2] = c->loc[2];
      temp->width[0] = c->width[0] / 2;
      temp->width[1] = c->width[1] / 2;
      temp->width[2] = c->width[2] / 2;
      temp->dmin = c->dmin / 2;
      if (k & 4) temp->loc[0] += temp->width[0];
      if (k & 2) temp->loc[1] += temp->width[1];
      if (k & 1) temp->loc[2] += temp->width[2];
      temp->depth = c->depth + 1;
      temp->split = 0;
      temp->hydro.dx_max_part = 0.f;
      temp->hydro.dx_max_sort = 0.f;
      temp->stars.dx_max_part = 0.f;
      temp->stars.dx_max_sort = 0.f;
      temp->nodeID = c->nodeID;
      temp->parent = c;
      c->progeny[k] = temp;
      c->split = 1;
      count += cell_unpack(&pc[pc->progeny[k]], temp, s, with_gravity);
    }

  /* Return the total number of unpacked cells. */
  c->mpi.pcell_size = count;
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Unpack the tags of a given cell and its sub-cells.
 *
 * @param tags An array of tags.
 * @param c The #cell in which to unpack the tags.
 *
 * @return The number of tags created.
 */
int cell_unpack_tags(const int *tags, struct cell *restrict c) {

#ifdef WITH_MPI

  /* Unpack the current pcell. */
  c->mpi.tag = tags[0];

  /* Number of new cells created. */
  int count = 1;

  /* Fill the progeny recursively, depth-first. */
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_unpack_tags(&tags[count], c->progeny[k]);
    }

#ifdef SWIFT_DEBUG_CHECKS
  if (c->mpi.pcell_size != count) error("Inconsistent tag and pcell count!");
#endif  // SWIFT_DEBUG_CHECKS

  /* Return the total number of unpacked tags. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Pack the time information of the given cell and all it's sub-cells.
 *
 * @param c The #cell.
 * @param pcells (output) The end-of-timestep information we pack into
 *
 * @return The number of packed cells.
 */
int cell_pack_end_step(struct cell *restrict c,
                       struct pcell_step *restrict pcells) {

#ifdef WITH_MPI

  /* Pack this cell's data. */
  pcells[0].hydro.ti_end_min = c->hydro.ti_end_min;
  pcells[0].hydro.ti_end_max = c->hydro.ti_end_max;
  pcells[0].grav.ti_end_min = c->grav.ti_end_min;
  pcells[0].grav.ti_end_max = c->grav.ti_end_max;
  pcells[0].stars.ti_end_min = c->stars.ti_end_min;
  pcells[0].stars.ti_end_max = c->stars.ti_end_max;
  pcells[0].hydro.dx_max_part = c->hydro.dx_max_part;
  pcells[0].stars.dx_max_part = c->stars.dx_max_part;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_pack_end_step(c->progeny[k], &pcells[count]);
    }

  /* Return the number of packed values. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Unpack the time information of a given cell and its sub-cells.
 *
 * @param c The #cell
 * @param pcells The end-of-timestep information to unpack
 *
 * @return The number of cells created.
 */
int cell_unpack_end_step(struct cell *restrict c,
                         struct pcell_step *restrict pcells) {

#ifdef WITH_MPI

  /* Unpack this cell's data. */
  c->hydro.ti_end_min = pcells[0].hydro.ti_end_min;
  c->hydro.ti_end_max = pcells[0].hydro.ti_end_max;
  c->grav.ti_end_min = pcells[0].grav.ti_end_min;
  c->grav.ti_end_max = pcells[0].grav.ti_end_max;
  c->stars.ti_end_min = pcells[0].stars.ti_end_min;
  c->stars.ti_end_max = pcells[0].stars.ti_end_max;
  c->hydro.dx_max_part = pcells[0].hydro.dx_max_part;
  c->stars.dx_max_part = pcells[0].stars.dx_max_part;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_unpack_end_step(c->progeny[k], &pcells[count]);
    }

  /* Return the number of packed values. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Pack the multipole information of the given cell and all it's
 * sub-cells.
 *
 * @param c The #cell.
 * @param pcells (output) The multipole information we pack into
 *
 * @return The number of packed cells.
 */
int cell_pack_multipoles(struct cell *restrict c,
                         struct gravity_tensors *restrict pcells) {

#ifdef WITH_MPI

  /* Pack this cell's data. */
  pcells[0] = *c->grav.multipole;

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_pack_multipoles(c->progeny[k], &pcells[count]);
    }

  /* Return the number of packed values. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Unpack the multipole information of a given cell and its sub-cells.
 *
 * @param c The #cell
 * @param pcells The multipole information to unpack
 *
 * @return The number of cells created.
 */
int cell_unpack_multipoles(struct cell *restrict c,
                           struct gravity_tensors *restrict pcells) {

#ifdef WITH_MPI

  /* Unpack this cell's data. */
  *c->grav.multipole = pcells[0];

  /* Fill in the progeny, depth-first recursion. */
  int count = 1;
  for (int k = 0; k < 8; k++)
    if (c->progeny[k] != NULL) {
      count += cell_unpack_multipoles(c->progeny[k], &pcells[count]);
    }

  /* Return the number of packed values. */
  return count;

#else
  error("SWIFT was not compiled with MPI support.");
  return 0;
#endif
}

/**
 * @brief Lock a cell for access to its array of #part and hold its parents.
 *
 * @param c The #cell.
 * @return 0 on success, 1 on failure
 */
int cell_locktree(struct cell *c) {

  TIMER_TIC

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

  TIMER_TIC

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

  TIMER_TIC

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

  TIMER_TIC

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
 * @brief Unlock a cell's parents for access to #part array.
 *
 * @param c The #cell.
 */
void cell_unlocktree(struct cell *c) {

  TIMER_TIC

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

  TIMER_TIC

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

  TIMER_TIC

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

  TIMER_TIC

  /* First of all, try to unlock this cell. */
  if (lock_unlock(&c->stars.lock) != 0) error("Failed to unlock cell.");

  /* Climb up the tree and unhold the parents. */
  for (struct cell *finger = c->parent; finger != NULL; finger = finger->parent)
    atomic_dec(&finger->stars.hold);

  TIMER_TOC(timer_locktree);
}

/**
 * @brief Sort the parts into eight bins along the given pivots.
 *
 * @param c The #cell array to be sorted.
 * @param parts_offset Offset of the cell parts array relative to the
 *        space's parts array, i.e. c->hydro.parts - s->parts.
 * @param sparts_offset Offset of the cell sparts array relative to the
 *        space's sparts array, i.e. c->stars.parts - s->stars.parts.
 * @param buff A buffer with at least max(c->hydro.count, c->grav.count)
 * entries, used for sorting indices.
 * @param sbuff A buffer with at least max(c->stars.count, c->grav.count)
 * entries, used for sorting indices for the sparts.
 * @param gbuff A buffer with at least max(c->hydro.count, c->grav.count)
 * entries, used for sorting indices for the gparts.
 */
void cell_split(struct cell *c, ptrdiff_t parts_offset, ptrdiff_t sparts_offset,
                struct cell_buff *buff, struct cell_buff *sbuff,
                struct cell_buff *gbuff) {

  const int count = c->hydro.count, gcount = c->grav.count,
            scount = c->stars.count;
  struct part *parts = c->hydro.parts;
  struct xpart *xparts = c->hydro.xparts;
  struct gpart *gparts = c->grav.parts;
  struct spart *sparts = c->stars.parts;
  const double pivot[3] = {c->loc[0] + c->width[0] / 2,
                           c->loc[1] + c->width[1] / 2,
                           c->loc[2] + c->width[2] / 2};
  int bucket_count[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int bucket_offset[9];

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that the buffs are OK. */
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
#endif /* SWIFT_DEBUG_CHECKS */

  /* Fill the buffer with the indices. */
  for (int k = 0; k < count; k++) {
    const int bid = (buff[k].x[0] >= pivot[0]) * 4 +
                    (buff[k].x[1] >= pivot[1]) * 2 + (buff[k].x[2] >= pivot[2]);
    bucket_count[bid]++;
    buff[k].ind = bid;
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
      int bid = buff[k].ind;
      if (bid != bucket) {
        struct part part = parts[k];
        struct xpart xpart = xparts[k];
        struct cell_buff temp_buff = buff[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (buff[j].ind == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&parts[j], &part, sizeof(struct part));
          memswap(&xparts[j], &xpart, sizeof(struct xpart));
          memswap(&buff[j], &temp_buff, sizeof(struct cell_buff));
          if (parts[j].gpart)
            parts[j].gpart->id_or_neg_offset = -(j + parts_offset);
          bid = temp_buff.ind;
        }
        parts[k] = part;
        xparts[k] = xpart;
        buff[k] = temp_buff;
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
    if (buff[k].ind < buff[k - 1].ind) error("Buff not sorted.");
    if (buff[k].x[0] != parts[k].x[0] || buff[k].x[1] != parts[k].x[1] ||
        buff[k].x[2] != parts[k].x[2])
      error("Inconsistent buff contents (k=%i).", k);
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
    const int bid = (sbuff[k].x[0] > pivot[0]) * 4 +
                    (sbuff[k].x[1] > pivot[1]) * 2 + (sbuff[k].x[2] > pivot[2]);
    bucket_count[bid]++;
    sbuff[k].ind = bid;
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
      int bid = sbuff[k].ind;
      if (bid != bucket) {
        struct spart spart = sparts[k];
        struct cell_buff temp_buff = sbuff[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (sbuff[j].ind == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&sparts[j], &spart, sizeof(struct spart));
          memswap(&sbuff[j], &temp_buff, sizeof(struct cell_buff));
          if (sparts[j].gpart)
            sparts[j].gpart->id_or_neg_offset = -(j + sparts_offset);
          bid = temp_buff.ind;
        }
        sparts[k] = spart;
        sbuff[k] = temp_buff;
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
  }

  /* Finally, do the same song and dance for the gparts. */
  for (int k = 0; k < 8; k++) bucket_count[k] = 0;

  /* Fill the buffer with the indices. */
  for (int k = 0; k < gcount; k++) {
    const int bid = (gbuff[k].x[0] > pivot[0]) * 4 +
                    (gbuff[k].x[1] > pivot[1]) * 2 + (gbuff[k].x[2] > pivot[2]);
    bucket_count[bid]++;
    gbuff[k].ind = bid;
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
      int bid = gbuff[k].ind;
      if (bid != bucket) {
        struct gpart gpart = gparts[k];
        struct cell_buff temp_buff = gbuff[k];
        while (bid != bucket) {
          int j = bucket_offset[bid] + bucket_count[bid]++;
          while (gbuff[j].ind == bid) {
            j++;
            bucket_count[bid]++;
          }
          memswap(&gparts[j], &gpart, sizeof(struct gpart));
          memswap(&gbuff[j], &temp_buff, sizeof(struct cell_buff));
          if (gparts[j].type == swift_type_gas) {
            parts[-gparts[j].id_or_neg_offset - parts_offset].gpart =
                &gparts[j];
          } else if (gparts[j].type == swift_type_stars) {
            sparts[-gparts[j].id_or_neg_offset - sparts_offset].gpart =
                &gparts[j];
          }
          bid = temp_buff.ind;
        }
        gparts[k] = gpart;
        gbuff[k] = temp_buff;
        if (gparts[k].type == swift_type_gas) {
          parts[-gparts[k].id_or_neg_offset - parts_offset].gpart = &gparts[k];
        } else if (gparts[k].type == swift_type_stars) {
          sparts[-gparts[k].id_or_neg_offset - sparts_offset].gpart =
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
  }
}

/**
 * @brief Sanitizes the smoothing length values of cells by setting large
 * outliers to more sensible values.
 *
 * Each cell with <1000 part will be processed. We limit h to be the size of
 * the cell and replace 0s with a good estimate.
 *
 * @param c The cell.
 * @param treated Has the cell already been sanitized at this level ?
 */
void cell_sanitize(struct cell *c, int treated) {

  const int count = c->hydro.count;
  const int scount = c->stars.count;
  struct part *parts = c->hydro.parts;
  struct spart *sparts = c->stars.parts;
  float h_max = 0.f;
  float stars_h_max = 0.f;

  /* Treat cells will <1000 particles */
  if (count < 1000 && !treated) {

    /* Get an upper bound on h */
    const float upper_h_max = c->dmin / (1.2f * kernel_gamma);

    /* Apply it */
    for (int i = 0; i < count; ++i) {
      if (parts[i].h == 0.f || parts[i].h > upper_h_max)
        parts[i].h = upper_h_max;
    }
    for (int i = 0; i < scount; ++i) {
      if (sparts[i].h == 0.f || sparts[i].h > upper_h_max)
        sparts[i].h = upper_h_max;
    }
  }

  /* Recurse and gather the new h_max values */
  if (c->split) {

    for (int k = 0; k < 8; ++k) {
      if (c->progeny[k] != NULL) {

        /* Recurse */
        cell_sanitize(c->progeny[k], (count < 1000));

        /* And collect */
        h_max = max(h_max, c->progeny[k]->hydro.h_max);
        stars_h_max = max(stars_h_max, c->progeny[k]->stars.h_max);
      }
    }
  } else {

    /* Get the new value of h_max */
    for (int i = 0; i < count; ++i) h_max = max(h_max, parts[i].h);
    for (int i = 0; i < scount; ++i)
      stars_h_max = max(stars_h_max, sparts[i].h);
  }

  /* Record the change */
  c->hydro.h_max = h_max;
  c->stars.h_max = stars_h_max;
}

/**
 * @brief Cleans the links in a given cell.
 *
 * @param c Cell to act upon
 * @param data Unused parameter
 */
void cell_clean_links(struct cell *c, void *data) {
  c->hydro.density = NULL;
  c->hydro.gradient = NULL;
  c->hydro.force = NULL;
  c->hydro.limiter = NULL;
  c->grav.grav = NULL;
  c->grav.mm = NULL;
  c->stars.density = NULL;
  c->stars.feedback = NULL;
}

/**
 * @brief Checks that the #part in a cell are at the
 * current point in time
 *
 * Calls error() if the cell is not at the current time.
 *
 * @param c Cell to act upon
 * @param data The current time on the integer time-line
 */
void cell_check_part_drift_point(struct cell *c, void *data) {

#ifdef SWIFT_DEBUG_CHECKS

  const integertime_t ti_drift = *(integertime_t *)data;

  /* Only check local cells */
  if (c->nodeID != engine_rank) return;

  /* Only check cells with content */
  if (c->hydro.count == 0) return;

  if (c->hydro.ti_old_part != ti_drift)
    error("Cell in an incorrect time-zone! c->hydro.ti_old=%lld ti_drift=%lld",
          c->hydro.ti_old_part, ti_drift);

  for (int i = 0; i < c->hydro.count; ++i)
    if (c->hydro.parts[i].ti_drift != ti_drift &&
        c->hydro.parts[i].time_bin != time_bin_inhibited)
      error("part in an incorrect time-zone! p->ti_drift=%lld ti_drift=%lld",
            c->hydro.parts[i].ti_drift, ti_drift);
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Checks that the #gpart in a cell are at the
 * current point in time
 *
 * Calls error() if the cell is not at the current time.
 *
 * @param c Cell to act upon
 * @param data The current time on the integer time-line
 */
void cell_check_gpart_drift_point(struct cell *c, void *data) {

#ifdef SWIFT_DEBUG_CHECKS

  const integertime_t ti_drift = *(integertime_t *)data;

  /* Only check local cells */
  if (c->nodeID != engine_rank) return;

  /* Only check cells with content */
  if (c->grav.count == 0) return;

  if (c->grav.ti_old_part != ti_drift)
    error(
        "Cell in an incorrect time-zone! c->grav.ti_old_part=%lld "
        "ti_drift=%lld",
        c->grav.ti_old_part, ti_drift);

  for (int i = 0; i < c->grav.count; ++i)
    if (c->grav.parts[i].ti_drift != ti_drift &&
        c->grav.parts[i].time_bin != time_bin_inhibited)
      error("g-part in an incorrect time-zone! gp->ti_drift=%lld ti_drift=%lld",
            c->grav.parts[i].ti_drift, ti_drift);
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Checks that the #spart in a cell are at the
 * current point in time
 *
 * Calls error() if the cell is not at the current time.
 *
 * @param c Cell to act upon
 * @param data The current time on the integer time-line
 */
void cell_check_spart_drift_point(struct cell *c, void *data) {

#ifdef SWIFT_DEBUG_CHECKS

  const integertime_t ti_drift = *(integertime_t *)data;

  /* Only check local cells */
  if (c->nodeID != engine_rank) return;

  /* Only check cells with content */
  if (c->stars.count == 0) return;

  if (c->stars.ti_old_part != ti_drift)
    error(
        "Cell in an incorrect time-zone! c->stars.ti_old_part=%lld "
        "ti_drift=%lld",
        c->stars.ti_old_part, ti_drift);

  for (int i = 0; i < c->stars.count; ++i)
    if (c->stars.parts[i].ti_drift != ti_drift &&
        c->stars.parts[i].time_bin != time_bin_inhibited)
      error("g-part in an incorrect time-zone! gp->ti_drift=%lld ti_drift=%lld",
            c->stars.parts[i].ti_drift, ti_drift);
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Checks that the multipole of a cell is at the current point in time
 *
 * Calls error() if the cell is not at the current time.
 *
 * @param c Cell to act upon
 * @param data The current time on the integer time-line
 */
void cell_check_multipole_drift_point(struct cell *c, void *data) {

#ifdef SWIFT_DEBUG_CHECKS

  const integertime_t ti_drift = *(integertime_t *)data;

  /* Only check local cells */
  if (c->nodeID != engine_rank) return;

  /* Only check cells with content */
  if (c->grav.count == 0) return;

  if (c->grav.ti_old_multipole != ti_drift)
    error(
        "Cell multipole in an incorrect time-zone! "
        "c->grav.ti_old_multipole=%lld "
        "ti_drift=%lld (depth=%d, node=%d)",
        c->grav.ti_old_multipole, ti_drift, c->depth, c->nodeID);

#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Resets all the individual cell task counters to 0.
 *
 * Should only be used for debugging purposes.
 *
 * @param c The #cell to reset.
 */
void cell_reset_task_counters(struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  for (int t = 0; t < task_type_count; ++t) c->tasks_executed[t] = 0;
  for (int t = 0; t < task_subtype_count; ++t) c->subtasks_executed[t] = 0;
  for (int k = 0; k < 8; ++k)
    if (c->progeny[k] != NULL) cell_reset_task_counters(c->progeny[k]);
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Recursively construct all the multipoles in a cell hierarchy.
 *
 * @param c The #cell.
 * @param ti_current The current integer time.
 */
void cell_make_multipoles(struct cell *c, integertime_t ti_current) {

  /* Reset everything */
  gravity_reset(c->grav.multipole);

  if (c->split) {

    /* Start by recursing */
    for (int k = 0; k < 8; ++k) {
      if (c->progeny[k] != NULL)
        cell_make_multipoles(c->progeny[k], ti_current);
    }

    /* Compute CoM of all progenies */
    double CoM[3] = {0., 0., 0.};
    double mass = 0.;

    for (int k = 0; k < 8; ++k) {
      if (c->progeny[k] != NULL) {
        const struct gravity_tensors *m = c->progeny[k]->grav.multipole;
        CoM[0] += m->CoM[0] * m->m_pole.M_000;
        CoM[1] += m->CoM[1] * m->m_pole.M_000;
        CoM[2] += m->CoM[2] * m->m_pole.M_000;
        mass += m->m_pole.M_000;
      }
    }

    const double mass_inv = 1. / mass;
    c->grav.multipole->CoM[0] = CoM[0] * mass_inv;
    c->grav.multipole->CoM[1] = CoM[1] * mass_inv;
    c->grav.multipole->CoM[2] = CoM[2] * mass_inv;

    /* Now shift progeny multipoles and add them up */
    struct multipole temp;
    double r_max = 0.;
    for (int k = 0; k < 8; ++k) {
      if (c->progeny[k] != NULL) {
        const struct cell *cp = c->progeny[k];
        const struct multipole *m = &cp->grav.multipole->m_pole;

        /* Contribution to multipole */
        gravity_M2M(&temp, m, c->grav.multipole->CoM, cp->grav.multipole->CoM);
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
    const double dx = c->grav.multipole->CoM[0] > c->loc[0] + c->width[0] * 0.5
                          ? c->grav.multipole->CoM[0] - c->loc[0]
                          : c->loc[0] + c->width[0] - c->grav.multipole->CoM[0];
    const double dy = c->grav.multipole->CoM[1] > c->loc[1] + c->width[1] * 0.5
                          ? c->grav.multipole->CoM[1] - c->loc[1]
                          : c->loc[1] + c->width[1] - c->grav.multipole->CoM[1];
    const double dz = c->grav.multipole->CoM[2] > c->loc[2] + c->width[2] * 0.5
                          ? c->grav.multipole->CoM[2] - c->loc[2]
                          : c->loc[2] + c->width[2] - c->grav.multipole->CoM[2];

    /* Take minimum of both limits */
    c->grav.multipole->r_max = min(r_max, sqrt(dx * dx + dy * dy + dz * dz));

  } else {

    if (c->grav.count > 0) {
      gravity_P2M(c->grav.multipole, c->grav.parts, c->grav.count);
      const double dx =
          c->grav.multipole->CoM[0] > c->loc[0] + c->width[0] * 0.5
              ? c->grav.multipole->CoM[0] - c->loc[0]
              : c->loc[0] + c->width[0] - c->grav.multipole->CoM[0];
      const double dy =
          c->grav.multipole->CoM[1] > c->loc[1] + c->width[1] * 0.5
              ? c->grav.multipole->CoM[1] - c->loc[1]
              : c->loc[1] + c->width[1] - c->grav.multipole->CoM[1];
      const double dz =
          c->grav.multipole->CoM[2] > c->loc[2] + c->width[2] * 0.5
              ? c->grav.multipole->CoM[2] - c->loc[2]
              : c->loc[2] + c->width[2] - c->grav.multipole->CoM[2];
      c->grav.multipole->r_max = sqrt(dx * dx + dy * dy + dz * dz);
    } else {
      gravity_multipole_init(&c->grav.multipole->m_pole);
      c->grav.multipole->CoM[0] = c->loc[0] + c->width[0] * 0.5;
      c->grav.multipole->CoM[1] = c->loc[1] + c->width[1] * 0.5;
      c->grav.multipole->CoM[2] = c->loc[2] + c->width[2] * 0.5;
      c->grav.multipole->r_max = 0.;
    }
  }

  /* Also update the values at rebuild time */
  c->grav.multipole->r_max_rebuild = c->grav.multipole->r_max;
  c->grav.multipole->CoM_rebuild[0] = c->grav.multipole->CoM[0];
  c->grav.multipole->CoM_rebuild[1] = c->grav.multipole->CoM[1];
  c->grav.multipole->CoM_rebuild[2] = c->grav.multipole->CoM[2];

  c->grav.ti_old_multipole = ti_current;
}

/**
 * @brief Recursively verify that the multipoles are the sum of their progenies.
 *
 * This function does not check whether the multipoles match the particle
 * content as we may not have received the particles.
 *
 * @param c The #cell to recursively search and verify.
 */
void cell_check_foreign_multipole(const struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS

  if (c->split) {

    double M_000 = 0.;
    long long num_gpart = 0;

    for (int k = 0; k < 8; k++) {
      const struct cell *cp = c->progeny[k];

      if (cp != NULL) {

        /* Check the mass */
        M_000 += cp->grav.multipole->m_pole.M_000;

        /* Check the number of particles */
        num_gpart += cp->grav.multipole->m_pole.num_gpart;

        /* Now recurse */
        cell_check_foreign_multipole(cp);
      }
    }

    if (num_gpart != c->grav.multipole->m_pole.num_gpart)
      error("Sum of particles in progenies does not match");
  }

#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Computes the multi-pole brutally and compare to the
 * recursively computed one.
 *
 * @param c Cell to act upon
 */
void cell_check_multipole(struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  struct gravity_tensors ma;
  const double tolerance = 1e-3; /* Relative */

  /* First recurse */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_check_multipole(c->progeny[k]);

  if (c->grav.count > 0) {

    /* Brute-force calculation */
    gravity_P2M(&ma, c->grav.parts, c->grav.count);

    /* Now  compare the multipole expansion */
    if (!gravity_multipole_equal(&ma, c->grav.multipole, tolerance)) {
      message("Multipoles are not equal at depth=%d! tol=%f", c->depth,
              tolerance);
      message("Correct answer:");
      gravity_multipole_print(&ma.m_pole);
      message("Recursive multipole:");
      gravity_multipole_print(&c->grav.multipole->m_pole);
      error("Aborting");
    }

    /* Check that the upper limit of r_max is good enough */
    if (!(1.1 * c->grav.multipole->r_max >= ma.r_max)) {
      error("Upper-limit r_max=%e too small. Should be >=%e.",
            c->grav.multipole->r_max, ma.r_max);
    } else if (c->grav.multipole->r_max * c->grav.multipole->r_max >
               3. * c->width[0] * c->width[0]) {
      error("r_max=%e larger than cell diagonal %e.", c->grav.multipole->r_max,
            sqrt(3. * c->width[0] * c->width[0]));
    }
  }
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Frees up the memory allocated for this #cell.
 *
 * @param c The #cell.
 */
void cell_clean(struct cell *c) {

  /* Hydro */
  cell_free_hydro_sorts(c);

  /* Stars */
  cell_free_stars_sorts(c);

  /* Recurse */
  for (int k = 0; k < 8; k++)
    if (c->progeny[k]) cell_clean(c->progeny[k]);
}

/**
 * @brief Clear the drift flags on the given cell.
 */
void cell_clear_drift_flags(struct cell *c, void *data) {
  c->hydro.do_drift = 0;
  c->hydro.do_sub_drift = 0;
  c->grav.do_drift = 0;
  c->grav.do_sub_drift = 0;
  c->stars.do_drift = 0;
  c->stars.do_sub_drift = 0;
}

/**
 * @brief Clear the limiter flags on the given cell.
 */
void cell_clear_limiter_flags(struct cell *c, void *data) {
  c->hydro.do_limiter = 0;
  c->hydro.do_sub_limiter = 0;
}

/**
 * @brief Activate the #part drifts on the given cell.
 */
void cell_activate_drift_part(struct cell *c, struct scheduler *s) {

  /* If this cell is already marked for drift, quit early. */
  if (c->hydro.do_drift) return;

  /* Mark this cell for drifting. */
  c->hydro.do_drift = 1;

  /* Set the do_sub_drifts all the way up and activate the super drift
     if this has not yet been done. */
  if (c == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->hydro.drift == NULL)
      error("Trying to activate un-existing c->hydro.drift");
#endif
    scheduler_activate(s, c->hydro.drift);
  } else {
    for (struct cell *parent = c->parent;
         parent != NULL && !parent->hydro.do_sub_drift;
         parent = parent->parent) {

      /* Mark this cell for drifting */
      parent->hydro.do_sub_drift = 1;

      if (parent == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
        if (parent->hydro.drift == NULL)
          error("Trying to activate un-existing parent->hydro.drift");
#endif
        scheduler_activate(s, parent->hydro.drift);
        break;
      }
    }
  }
}

/**
 * @brief Activate the #gpart drifts on the given cell.
 */
void cell_activate_drift_gpart(struct cell *c, struct scheduler *s) {

  /* If this cell is already marked for drift, quit early. */
  if (c->grav.do_drift) return;

  /* Mark this cell for drifting. */
  c->grav.do_drift = 1;

  if (c->grav.drift_out != NULL) scheduler_activate(s, c->grav.drift_out);

  /* Set the do_grav_sub_drifts all the way up and activate the super drift
     if this has not yet been done. */
  if (c == c->grav.super) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->grav.drift == NULL)
      error("Trying to activate un-existing c->grav.drift");
#endif
    scheduler_activate(s, c->grav.drift);
  } else {
    for (struct cell *parent = c->parent;
         parent != NULL && !parent->grav.do_sub_drift;
         parent = parent->parent) {
      parent->grav.do_sub_drift = 1;

      if (parent->grav.drift_out) {
        scheduler_activate(s, parent->grav.drift_out);
      }

      if (parent == c->grav.super) {
#ifdef SWIFT_DEBUG_CHECKS
        if (parent->grav.drift == NULL)
          error("Trying to activate un-existing parent->grav.drift");
#endif
        scheduler_activate(s, parent->grav.drift);
        break;
      }
    }
  }
}

/**
 * @brief Activate the #spart drifts on the given cell.
 */
void cell_activate_drift_spart(struct cell *c, struct scheduler *s) {

  /* If this cell is already marked for drift, quit early. */
  if (c->stars.do_drift) return;

  /* Mark this cell for drifting. */
  c->stars.do_drift = 1;

  /* Set the do_stars_sub_drifts all the way up and activate the super drift
     if this has not yet been done. */
  if (c == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->stars.drift == NULL)
      error("Trying to activate un-existing c->stars.drift");
#endif
    scheduler_activate(s, c->stars.drift);
  } else {
    for (struct cell *parent = c->parent;
         parent != NULL && !parent->stars.do_sub_drift;
         parent = parent->parent) {

      /* Mark this cell for drifting */
      parent->stars.do_sub_drift = 1;

      if (parent == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
        if (parent->stars.drift == NULL)
          error("Trying to activate un-existing parent->stars.drift");
#endif
        scheduler_activate(s, parent->stars.drift);
        break;
      }
    }
  }
}

/**
 * @brief Activate the drifts on the given cell.
 */
void cell_activate_limiter(struct cell *c, struct scheduler *s) {

  /* If this cell is already marked for limiting, quit early. */
  if (c->hydro.do_limiter) return;

  /* Mark this cell for limiting. */
  c->hydro.do_limiter = 1;

  /* Set the do_sub_limiter all the way up and activate the super limiter
     if this has not yet been done. */
  if (c == c->super) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->timestep_limiter == NULL)
      error("Trying to activate un-existing c->timestep_limiter");
#endif
    scheduler_activate(s, c->timestep_limiter);
  } else {
    for (struct cell *parent = c->parent;
         parent != NULL && !parent->hydro.do_sub_limiter;
         parent = parent->parent) {

      /* Mark this cell for limiting */
      parent->hydro.do_sub_limiter = 1;

      if (parent == c->super) {
#ifdef SWIFT_DEBUG_CHECKS
        if (parent->timestep_limiter == NULL)
          error("Trying to activate un-existing parent->timestep_limiter");
#endif
        scheduler_activate(s, parent->timestep_limiter);
        break;
      }
    }
  }
}

/**
 * @brief Activate the sorts up a cell hierarchy.
 */
void cell_activate_hydro_sorts_up(struct cell *c, struct scheduler *s) {

  if (c == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->hydro.sorts == NULL)
      error("Trying to activate un-existing c->hydro.sorts");
#endif
    scheduler_activate(s, c->hydro.sorts);
    if (c->nodeID == engine_rank) cell_activate_drift_part(c, s);
  } else {

    for (struct cell *parent = c->parent;
         parent != NULL && !parent->hydro.do_sub_sort;
         parent = parent->parent) {
      parent->hydro.do_sub_sort = 1;
      if (parent == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
        if (parent->hydro.sorts == NULL)
          error("Trying to activate un-existing parents->hydro.sorts");
#endif
        scheduler_activate(s, parent->hydro.sorts);
        if (parent->nodeID == engine_rank) cell_activate_drift_part(parent, s);
        break;
      }
    }
  }
}

/**
 * @brief Activate the sorts on a given cell, if needed.
 */
void cell_activate_hydro_sorts(struct cell *c, int sid, struct scheduler *s) {

  /* Do we need to re-sort? */
  if (c->hydro.dx_max_sort > space_maxreldx * c->dmin) {

    /* Climb up the tree to active the sorts in that direction */
    for (struct cell *finger = c; finger != NULL; finger = finger->parent) {
      if (finger->hydro.requires_sorts) {
        atomic_or(&finger->hydro.do_sort, finger->hydro.requires_sorts);
        cell_activate_hydro_sorts_up(finger, s);
      }
      finger->hydro.sorted = 0;
    }
  }

  /* Has this cell been sorted at all for the given sid? */
  if (!(c->hydro.sorted & (1 << sid)) || c->nodeID != engine_rank) {

    atomic_or(&c->hydro.do_sort, (1 << sid));
    cell_activate_hydro_sorts_up(c, s);
  }
}

/**
 * @brief Activate the sorts up a cell hierarchy.
 */
void cell_activate_stars_sorts_up(struct cell *c, struct scheduler *s) {

  if (c == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->stars.sorts == NULL)
      error("Trying to activate un-existing c->stars.sorts");
#endif
    scheduler_activate(s, c->stars.sorts);
    if (c->nodeID == engine_rank) {
      // MATTHIEU: to do: do we actually need both drifts here?
      cell_activate_drift_spart(c, s);
    }
  } else {

    for (struct cell *parent = c->parent;
         parent != NULL && !parent->stars.do_sub_sort;
         parent = parent->parent) {
      parent->stars.do_sub_sort = 1;
      if (parent == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
        if (parent->stars.sorts == NULL)
          error("Trying to activate un-existing parents->stars.sorts");
#endif
        scheduler_activate(s, parent->stars.sorts);
        if (parent->nodeID == engine_rank) cell_activate_drift_spart(parent, s);
        break;
      }
    }
  }
}

/**
 * @brief Activate the sorts on a given cell, if needed.
 */
void cell_activate_stars_sorts(struct cell *c, int sid, struct scheduler *s) {

  /* Do we need to re-sort? */
  if (c->stars.dx_max_sort > space_maxreldx * c->dmin) {

    /* Climb up the tree to active the sorts in that direction */
    for (struct cell *finger = c; finger != NULL; finger = finger->parent) {
      if (finger->stars.requires_sorts) {
        atomic_or(&finger->stars.do_sort, finger->stars.requires_sorts);
        cell_activate_stars_sorts_up(finger, s);
      }
      finger->stars.sorted = 0;
    }
  }

  /* Has this cell been sorted at all for the given sid? */
  if (!(c->stars.sorted & (1 << sid)) || c->nodeID != engine_rank) {

    atomic_or(&c->stars.do_sort, (1 << sid));
    cell_activate_stars_sorts_up(c, s);
  }
}

/**
 * @brief Traverse a sub-cell task and activate the hydro drift tasks that are
 * required by a hydro task
 *
 * @param ci The first #cell we recurse in.
 * @param cj The second #cell we recurse in.
 * @param s The task #scheduler.
 */
void cell_activate_subcell_hydro_tasks(struct cell *ci, struct cell *cj,
                                       struct scheduler *s) {
  const struct engine *e = s->space->e;
  const int with_limiter = (e->policy & engine_policy_limiter);

  /* Store the current dx_max and h_max values. */
  ci->hydro.dx_max_part_old = ci->hydro.dx_max_part;
  ci->hydro.h_max_old = ci->hydro.h_max;

  if (cj != NULL) {
    cj->hydro.dx_max_part_old = cj->hydro.dx_max_part;
    cj->hydro.h_max_old = cj->hydro.h_max;
  }

  /* Self interaction? */
  if (cj == NULL) {

    /* Do anything? */
    if (ci->hydro.count == 0 || !cell_is_active_hydro(ci, e)) return;

    /* Recurse? */
    if (cell_can_recurse_in_self_hydro_task(ci)) {

      /* Loop over all progenies and pairs of progenies */
      for (int j = 0; j < 8; j++) {
        if (ci->progeny[j] != NULL) {
          cell_activate_subcell_hydro_tasks(ci->progeny[j], NULL, s);
          for (int k = j + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              cell_activate_subcell_hydro_tasks(ci->progeny[j], ci->progeny[k],
                                                s);
        }
      }
    } else {

      /* We have reached the bottom of the tree: activate drift */
      cell_activate_drift_part(ci, s);
      if (with_limiter) cell_activate_limiter(ci, s);
    }
  }

  /* Otherwise, pair interation */
  else {

    /* Should we even bother? */
    if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;
    if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

    /* Get the orientation of the pair. */
    double shift[3];
    int sid = space_getsid(s->space, &ci, &cj, shift);

    /* recurse? */
    if (cell_can_recurse_in_pair_hydro_task(ci) &&
        cell_can_recurse_in_pair_hydro_task(cj)) {

      /* Different types of flags. */
      switch (sid) {

        /* Regular sub-cell interactions of a single cell. */
        case 0: /* (  1 ,  1 ,  1 ) */
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          break;

        case 1: /* (  1 ,  1 ,  0 ) */
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[0],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[1],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[1],
                                              s);
          break;

        case 2: /* (  1 ,  1 , -1 ) */
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[1],
                                              s);
          break;

        case 3: /* (  1 ,  0 ,  1 ) */
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[0],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[2],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[2],
                                              s);
          break;

        case 4: /* (  1 ,  0 ,  0 ) */
          if (ci->progeny[4] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[4], cj->progeny[0],
                                              s);
          if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[4], cj->progeny[1],
                                              s);
          if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[4], cj->progeny[2],
                                              s);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[4], cj->progeny[3],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[0],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[1],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[2],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[3],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[0],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[1],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[2],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[3],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[1],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[2],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[3],
                                              s);
          break;

        case 5: /* (  1 ,  0 , -1 ) */
          if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[4], cj->progeny[1],
                                              s);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[4], cj->progeny[3],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[1],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[3],
                                              s);
          break;

        case 6: /* (  1 , -1 ,  1 ) */
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[2],
                                              s);
          break;

        case 7: /* (  1 , -1 ,  0 ) */
          if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[4], cj->progeny[2],
                                              s);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[4], cj->progeny[3],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[2],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[3],
                                              s);
          break;

        case 8: /* (  1 , -1 , -1 ) */
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[4], cj->progeny[3],
                                              s);
          break;

        case 9: /* (  0 ,  1 ,  1 ) */
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[3], cj->progeny[0],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[3], cj->progeny[4],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[4],
                                              s);
          break;

        case 10: /* (  0 ,  1 ,  0 ) */
          if (ci->progeny[2] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[2], cj->progeny[0],
                                              s);
          if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[2], cj->progeny[1],
                                              s);
          if (ci->progeny[2] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[2], cj->progeny[4],
                                              s);
          if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[2], cj->progeny[5],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[3], cj->progeny[0],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[3], cj->progeny[1],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[3], cj->progeny[4],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[5] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[3], cj->progeny[5],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[0],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[1],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[4],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[5],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[1],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[4],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[5] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[5],
                                              s);
          break;

        case 11: /* (  0 ,  1 , -1 ) */
          if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[2], cj->progeny[1],
                                              s);
          if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[2], cj->progeny[5],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[1],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[6], cj->progeny[5],
                                              s);
          break;

        case 12: /* (  0 ,  0 ,  1 ) */
          if (ci->progeny[1] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[1], cj->progeny[0],
                                              s);
          if (ci->progeny[1] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[1], cj->progeny[2],
                                              s);
          if (ci->progeny[1] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[1], cj->progeny[4],
                                              s);
          if (ci->progeny[1] != NULL && cj->progeny[6] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[1], cj->progeny[6],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[3], cj->progeny[0],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[3], cj->progeny[2],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[3], cj->progeny[4],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[6] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[3], cj->progeny[6],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[0],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[2],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[4],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[6] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[5], cj->progeny[6],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[2],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[4],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[6] != NULL)
            cell_activate_subcell_hydro_tasks(ci->progeny[7], cj->progeny[6],
                                              s);
          break;
      }

    }

    /* Otherwise, activate the sorts and drifts. */
    else if (cell_is_active_hydro(ci, e) || cell_is_active_hydro(cj, e)) {

      /* We are going to interact this pair, so store some values. */
      atomic_or(&ci->hydro.requires_sorts, 1 << sid);
      atomic_or(&cj->hydro.requires_sorts, 1 << sid);
      ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;
      cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;

      /* Activate the drifts if the cells are local. */
      if (ci->nodeID == engine_rank) cell_activate_drift_part(ci, s);
      if (cj->nodeID == engine_rank) cell_activate_drift_part(cj, s);

      /* Also activate the time-step limiter */
      if (ci->nodeID == engine_rank && with_limiter)
        cell_activate_limiter(ci, s);
      if (cj->nodeID == engine_rank && with_limiter)
        cell_activate_limiter(cj, s);

      /* Do we need to sort the cells? */
      cell_activate_hydro_sorts(ci, sid, s);
      cell_activate_hydro_sorts(cj, sid, s);
    }
  } /* Otherwise, pair interation */
}

/**
 * @brief Traverse a sub-cell task and activate the stars drift tasks that are
 * required by a stars task
 *
 * @param ci The first #cell we recurse in.
 * @param cj The second #cell we recurse in.
 * @param s The task #scheduler.
 */
void cell_activate_subcell_stars_tasks(struct cell *ci, struct cell *cj,
                                       struct scheduler *s) {
  const struct engine *e = s->space->e;

  /* Store the current dx_max and h_max values. */
  ci->stars.dx_max_part_old = ci->stars.dx_max_part;
  ci->stars.h_max_old = ci->stars.h_max;
  ci->hydro.dx_max_part_old = ci->hydro.dx_max_part;
  ci->hydro.h_max_old = ci->hydro.h_max;

  if (cj != NULL) {
    cj->stars.dx_max_part_old = cj->stars.dx_max_part;
    cj->stars.h_max_old = cj->stars.h_max;
    cj->hydro.dx_max_part_old = cj->hydro.dx_max_part;
    cj->hydro.h_max_old = cj->hydro.h_max;
  }

  /* Self interaction? */
  if (cj == NULL) {

    /* Do anything? */
    if (!cell_is_active_stars(ci, e) || ci->hydro.count == 0 ||
        ci->stars.count == 0)
      return;

    /* Recurse? */
    if (cell_can_recurse_in_self_stars_task(ci)) {

      /* Loop over all progenies and pairs of progenies */
      for (int j = 0; j < 8; j++) {
        if (ci->progeny[j] != NULL) {
          cell_activate_subcell_stars_tasks(ci->progeny[j], NULL, s);
          for (int k = j + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              cell_activate_subcell_stars_tasks(ci->progeny[j], ci->progeny[k],
                                                s);
        }
      }
    } else {

      /* We have reached the bottom of the tree: activate drift */
      cell_activate_drift_spart(ci, s);
      cell_activate_drift_part(ci, s);
    }
  }

  /* Otherwise, pair interation */
  else {

    /* Should we even bother? */
    if (!cell_is_active_stars(ci, e) && !cell_is_active_stars(cj, e)) return;

    /* Get the orientation of the pair. */
    double shift[3];
    int sid = space_getsid(s->space, &ci, &cj, shift);

    /* recurse? */
    if (cell_can_recurse_in_pair_stars_task(ci, cj) &&
        cell_can_recurse_in_pair_stars_task(cj, ci)) {

      /* Different types of flags. */
      switch (sid) {

        /* Regular sub-cell interactions of a single cell. */
        case 0: /* (  1 ,  1 ,  1 ) */
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          break;

        case 1: /* (  1 ,  1 ,  0 ) */
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[0],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[1],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[1],
                                              s);
          break;

        case 2: /* (  1 ,  1 , -1 ) */
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[1],
                                              s);
          break;

        case 3: /* (  1 ,  0 ,  1 ) */
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[0],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[2],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[2],
                                              s);
          break;

        case 4: /* (  1 ,  0 ,  0 ) */
          if (ci->progeny[4] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[4], cj->progeny[0],
                                              s);
          if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[4], cj->progeny[1],
                                              s);
          if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[4], cj->progeny[2],
                                              s);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[4], cj->progeny[3],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[0],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[1],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[2],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[3],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[0],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[1],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[2],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[3],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[1],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[2],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[3],
                                              s);
          break;

        case 5: /* (  1 ,  0 , -1 ) */
          if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[4], cj->progeny[1],
                                              s);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[4], cj->progeny[3],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[1],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[3],
                                              s);
          break;

        case 6: /* (  1 , -1 ,  1 ) */
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[2],
                                              s);
          break;

        case 7: /* (  1 , -1 ,  0 ) */
          if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[4], cj->progeny[2],
                                              s);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[4], cj->progeny[3],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[2],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[3],
                                              s);
          break;

        case 8: /* (  1 , -1 , -1 ) */
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[4], cj->progeny[3],
                                              s);
          break;

        case 9: /* (  0 ,  1 ,  1 ) */
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[3], cj->progeny[0],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[3], cj->progeny[4],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[4],
                                              s);
          break;

        case 10: /* (  0 ,  1 ,  0 ) */
          if (ci->progeny[2] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[2], cj->progeny[0],
                                              s);
          if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[2], cj->progeny[1],
                                              s);
          if (ci->progeny[2] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[2], cj->progeny[4],
                                              s);
          if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[2], cj->progeny[5],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[3], cj->progeny[0],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[3], cj->progeny[1],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[3], cj->progeny[4],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[5] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[3], cj->progeny[5],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[0],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[1],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[4],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[5],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[1],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[4],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[5] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[5],
                                              s);
          break;

        case 11: /* (  0 ,  1 , -1 ) */
          if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[2], cj->progeny[1],
                                              s);
          if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[2], cj->progeny[5],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[1],
                                              s);
          if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[6], cj->progeny[5],
                                              s);
          break;

        case 12: /* (  0 ,  0 ,  1 ) */
          if (ci->progeny[1] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[1], cj->progeny[0],
                                              s);
          if (ci->progeny[1] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[1], cj->progeny[2],
                                              s);
          if (ci->progeny[1] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[1], cj->progeny[4],
                                              s);
          if (ci->progeny[1] != NULL && cj->progeny[6] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[1], cj->progeny[6],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[3], cj->progeny[0],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[3], cj->progeny[2],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[3], cj->progeny[4],
                                              s);
          if (ci->progeny[3] != NULL && cj->progeny[6] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[3], cj->progeny[6],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[0],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[2],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[4],
                                              s);
          if (ci->progeny[5] != NULL && cj->progeny[6] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[5], cj->progeny[6],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[0],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[2],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[4],
                                              s);
          if (ci->progeny[7] != NULL && cj->progeny[6] != NULL)
            cell_activate_subcell_stars_tasks(ci->progeny[7], cj->progeny[6],
                                              s);
          break;
      }

    }

    /* Otherwise, activate the sorts and drifts. */
    else {

      if (cell_is_active_stars(ci, e)) {
        /* We are going to interact this pair, so store some values. */
        atomic_or(&cj->hydro.requires_sorts, 1 << sid);
        atomic_or(&ci->stars.requires_sorts, 1 << sid);

        cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;
        ci->stars.dx_max_sort_old = ci->stars.dx_max_sort;

        /* Activate the drifts if the cells are local. */
        if (ci->nodeID == engine_rank) cell_activate_drift_spart(ci, s);
        if (cj->nodeID == engine_rank) cell_activate_drift_part(cj, s);

        /* Do we need to sort the cells? */
        cell_activate_hydro_sorts(cj, sid, s);
        cell_activate_stars_sorts(ci, sid, s);
      }

      if (cell_is_active_stars(cj, e)) {
        /* We are going to interact this pair, so store some values. */
        atomic_or(&cj->stars.requires_sorts, 1 << sid);
        atomic_or(&ci->hydro.requires_sorts, 1 << sid);

        ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;
        cj->stars.dx_max_sort_old = cj->stars.dx_max_sort;

        /* Activate the drifts if the cells are local. */
        if (ci->nodeID == engine_rank) cell_activate_drift_part(ci, s);
        if (cj->nodeID == engine_rank) cell_activate_drift_spart(cj, s);

        /* Do we need to sort the cells? */
        cell_activate_hydro_sorts(ci, sid, s);
        cell_activate_stars_sorts(cj, sid, s);
      }
    }
  } /* Otherwise, pair interation */
}

/**
 * @brief Traverse a sub-cell task and activate the gravity drift tasks that
 * are required by a self gravity task.
 *
 * @param ci The first #cell we recurse in.
 * @param cj The second #cell we recurse in.
 * @param s The task #scheduler.
 */
void cell_activate_subcell_grav_tasks(struct cell *ci, struct cell *cj,
                                      struct scheduler *s) {
  /* Some constants */
  const struct space *sp = s->space;
  const struct engine *e = sp->e;

  /* Self interaction? */
  if (cj == NULL) {

    /* Do anything? */
    if (ci->grav.count == 0 || !cell_is_active_gravity(ci, e)) return;

    /* Recurse? */
    if (ci->split) {

      /* Loop over all progenies and pairs of progenies */
      for (int j = 0; j < 8; j++) {
        if (ci->progeny[j] != NULL) {
          cell_activate_subcell_grav_tasks(ci->progeny[j], NULL, s);
          for (int k = j + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              cell_activate_subcell_grav_tasks(ci->progeny[j], ci->progeny[k],
                                               s);
        }
      }
    } else {

      /* We have reached the bottom of the tree: activate gpart drift */
      cell_activate_drift_gpart(ci, s);
    }
  }

  /* Pair interaction */
  else {

    /* Anything to do here? */
    if (!cell_is_active_gravity(ci, e) && !cell_is_active_gravity(cj, e))
      return;
    if (ci->grav.count == 0 || cj->grav.count == 0) return;

    /* Atomically drift the multipole in ci */
    lock_lock(&ci->grav.mlock);
    if (ci->grav.ti_old_multipole < e->ti_current) cell_drift_multipole(ci, e);
    if (lock_unlock(&ci->grav.mlock) != 0) error("Impossible to unlock m-pole");

    /* Atomically drift the multipole in cj */
    lock_lock(&cj->grav.mlock);
    if (cj->grav.ti_old_multipole < e->ti_current) cell_drift_multipole(cj, e);
    if (lock_unlock(&cj->grav.mlock) != 0) error("Impossible to unlock m-pole");

    /* Can we use multipoles ? */
    if (cell_can_use_pair_mm(ci, cj, e, sp)) {

      /* Ok, no need to drift anything */
      return;
    }
    /* Otherwise, activate the gpart drifts if we are at the bottom. */
    else if (!ci->split && !cj->split) {

      /* Activate the drifts if the cells are local. */
      if (cell_is_active_gravity(ci, e) || cell_is_active_gravity(cj, e)) {
        if (ci->nodeID == engine_rank) cell_activate_drift_gpart(ci, s);
        if (cj->nodeID == engine_rank) cell_activate_drift_gpart(cj, s);
      }
    }
    /* Ok, we can still recurse */
    else {

      /* Recover the multipole information */
      const struct gravity_tensors *const multi_i = ci->grav.multipole;
      const struct gravity_tensors *const multi_j = cj->grav.multipole;
      const double ri_max = multi_i->r_max;
      const double rj_max = multi_j->r_max;

      if (ri_max > rj_max) {
        if (ci->split) {

          /* Loop over ci's children */
          for (int k = 0; k < 8; k++) {
            if (ci->progeny[k] != NULL)
              cell_activate_subcell_grav_tasks(ci->progeny[k], cj, s);
          }

        } else if (cj->split) {

          /* Loop over cj's children */
          for (int k = 0; k < 8; k++) {
            if (cj->progeny[k] != NULL)
              cell_activate_subcell_grav_tasks(ci, cj->progeny[k], s);
          }

        } else {
          error("Fundamental error in the logic");
        }
      } else if (rj_max >= ri_max) {
        if (cj->split) {

          /* Loop over cj's children */
          for (int k = 0; k < 8; k++) {
            if (cj->progeny[k] != NULL)
              cell_activate_subcell_grav_tasks(ci, cj->progeny[k], s);
          }

        } else if (ci->split) {

          /* Loop over ci's children */
          for (int k = 0; k < 8; k++) {
            if (ci->progeny[k] != NULL)
              cell_activate_subcell_grav_tasks(ci->progeny[k], cj, s);
          }

        } else {
          error("Fundamental error in the logic");
        }
      }
    }
  }
}

/**
 * @brief Traverse a sub-cell task and activate the gravity drift tasks that
 * are required by an external gravity task.
 *
 * @param ci The #cell we recurse in.
 * @param s The task #scheduler.
 */
void cell_activate_subcell_external_grav_tasks(struct cell *ci,
                                               struct scheduler *s) {

  /* Some constants */
  const struct space *sp = s->space;
  const struct engine *e = sp->e;

  /* Do anything? */
  if (!cell_is_active_gravity(ci, e)) return;

  /* Recurse? */
  if (ci->split) {

    /* Loop over all progenies (no need for pairs for self-gravity) */
    for (int j = 0; j < 8; j++) {
      if (ci->progeny[j] != NULL) {
        cell_activate_subcell_external_grav_tasks(ci->progeny[j], s);
      }
    }
  } else {

    /* We have reached the bottom of the tree: activate gpart drift */
    cell_activate_drift_gpart(ci, s);
  }
}

/**
 * @brief Un-skips all the hydro tasks associated with a given cell and checks
 * if the space needs to be rebuilt.
 *
 * @param c the #cell.
 * @param s the #scheduler.
 *
 * @return 1 If the space needs rebuilding. 0 otherwise.
 */
int cell_unskip_hydro_tasks(struct cell *c, struct scheduler *s) {

  struct engine *e = s->space->e;
  const int nodeID = e->nodeID;
  const int with_limiter = (e->policy & engine_policy_limiter);
  int rebuild = 0;

  /* Un-skip the density tasks involved with this cell. */
  for (struct link *l = c->hydro.density; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const int ci_active = cell_is_active_hydro(ci, e);
    const int cj_active = (cj != NULL) ? cell_is_active_hydro(cj, e) : 0;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active && ci_nodeID == nodeID) ||
        (cj_active && cj_nodeID == nodeID)) {
      scheduler_activate(s, t);

      /* Activate hydro drift */
      if (t->type == task_type_self) {
        if (ci_nodeID == nodeID) cell_activate_drift_part(ci, s);
        if (ci_nodeID == nodeID && with_limiter) cell_activate_limiter(ci, s);
      }

      /* Set the correct sorting flags and activate hydro drifts */
      else if (t->type == task_type_pair) {
        /* Store some values. */
        atomic_or(&ci->hydro.requires_sorts, 1 << t->flags);
        atomic_or(&cj->hydro.requires_sorts, 1 << t->flags);
        ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;
        cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;

        /* Activate the drift tasks. */
        if (ci_nodeID == nodeID) cell_activate_drift_part(ci, s);
        if (cj_nodeID == nodeID) cell_activate_drift_part(cj, s);

        /* Activate the limiter tasks. */
        if (ci_nodeID == nodeID && with_limiter) cell_activate_limiter(ci, s);
        if (cj_nodeID == nodeID && with_limiter) cell_activate_limiter(cj, s);

        /* Check the sorts and activate them if needed. */
        cell_activate_hydro_sorts(ci, t->flags, s);
        cell_activate_hydro_sorts(cj, t->flags, s);
      }
      /* Store current values of dx_max and h_max. */
      else if (t->type == task_type_sub_pair || t->type == task_type_sub_self) {
        cell_activate_subcell_hydro_tasks(ci, cj, s);
      }
    }

    /* Only interested in pair interactions as of here. */
    if (t->type == task_type_pair || t->type == task_type_sub_pair) {

      /* Check whether there was too much particle motion, i.e. the
         cell neighbour conditions were violated. */
      if (cell_need_rebuild_for_hydro_pair(ci, cj)) rebuild = 1;

#ifdef WITH_MPI
      /* Activate the send/recv tasks. */
      if (ci_nodeID != nodeID) {

        /* If the local cell is active, receive data from the foreign cell. */
        if (cj_active) {
          scheduler_activate(s, ci->mpi.hydro.recv_xv);
          if (ci_active) {
            scheduler_activate(s, ci->mpi.hydro.recv_rho);

#ifdef EXTRA_HYDRO_LOOP
            scheduler_activate(s, ci->mpi.hydro.recv_gradient);
#endif
          }
        }

        /* If the foreign cell is active, we want its ti_end values. */
        if (ci_active || with_limiter) scheduler_activate(s, ci->mpi.recv_ti);

        if (with_limiter) scheduler_activate(s, ci->mpi.limiter.recv);
        if (with_limiter)
          scheduler_activate_send(s, cj->mpi.limiter.send, ci->nodeID);

        /* Is the foreign cell active and will need stuff from us? */
        if (ci_active) {

          scheduler_activate_send(s, cj->mpi.hydro.send_xv, ci_nodeID);

          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(cj, s);
          if (with_limiter) cell_activate_limiter(cj, s);

          /* If the local cell is also active, more stuff will be needed. */
          if (cj_active) {
            scheduler_activate_send(s, cj->mpi.hydro.send_rho, ci_nodeID);

#ifdef EXTRA_HYDRO_LOOP
            scheduler_activate_send(s, cj->mpi.hydro.send_gradient, ci_nodeID);
#endif
          }
        }

        /* If the local cell is active, send its ti_end values. */
        if (cj_active || with_limiter)
          scheduler_activate_send(s, cj->mpi.send_ti, ci_nodeID);

      } else if (cj_nodeID != nodeID) {

        /* If the local cell is active, receive data from the foreign cell. */
        if (ci_active) {
          scheduler_activate(s, cj->mpi.hydro.recv_xv);
          if (cj_active) {
            scheduler_activate(s, cj->mpi.hydro.recv_rho);

#ifdef EXTRA_HYDRO_LOOP
            scheduler_activate(s, cj->mpi.hydro.recv_gradient);
#endif
          }
        }

        /* If the foreign cell is active, we want its ti_end values. */
        if (cj_active || with_limiter) scheduler_activate(s, cj->mpi.recv_ti);

        if (with_limiter) scheduler_activate(s, cj->mpi.limiter.recv);
        if (with_limiter)
          scheduler_activate_send(s, ci->mpi.limiter.send, cj->nodeID);

        /* Is the foreign cell active and will need stuff from us? */
        if (cj_active) {

          scheduler_activate_send(s, ci->mpi.hydro.send_xv, cj_nodeID);

          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(ci, s);
          if (with_limiter) cell_activate_limiter(ci, s);

          /* If the local cell is also active, more stuff will be needed. */
          if (ci_active) {

            scheduler_activate_send(s, ci->mpi.hydro.send_rho, cj_nodeID);

#ifdef EXTRA_HYDRO_LOOP
            scheduler_activate_send(s, ci->mpi.hydro.send_gradient, cj_nodeID);
#endif
          }
        }

        /* If the local cell is active, send its ti_end values. */
        if (ci_active || with_limiter)
          scheduler_activate_send(s, ci->mpi.send_ti, cj_nodeID);
      }
#endif
    }
  }

  /* Unskip all the other task types. */
  if (c->nodeID == nodeID && cell_is_active_hydro(c, e)) {

    for (struct link *l = c->hydro.gradient; l != NULL; l = l->next)
      scheduler_activate(s, l->t);
    for (struct link *l = c->hydro.force; l != NULL; l = l->next)
      scheduler_activate(s, l->t);
    for (struct link *l = c->hydro.limiter; l != NULL; l = l->next)
      scheduler_activate(s, l->t);

    if (c->hydro.extra_ghost != NULL)
      scheduler_activate(s, c->hydro.extra_ghost);
    if (c->hydro.ghost_in != NULL) scheduler_activate(s, c->hydro.ghost_in);
    if (c->hydro.ghost_out != NULL) scheduler_activate(s, c->hydro.ghost_out);
    if (c->hydro.ghost != NULL) scheduler_activate(s, c->hydro.ghost);
    if (c->kick1 != NULL) scheduler_activate(s, c->kick1);
    if (c->kick2 != NULL) scheduler_activate(s, c->kick2);
    if (c->timestep != NULL) scheduler_activate(s, c->timestep);
    if (c->hydro.end_force != NULL) scheduler_activate(s, c->hydro.end_force);
    if (c->hydro.cooling != NULL) scheduler_activate(s, c->hydro.cooling);
    if (c->logger != NULL) scheduler_activate(s, c->logger);

    if (c->hydro.star_formation != NULL) {
      scheduler_activate(s, c->hydro.star_formation);
      cell_activate_drift_spart(c, s);
    }
  }

  return rebuild;
}

/**
 * @brief Un-skips all the gravity tasks associated with a given cell and checks
 * if the space needs to be rebuilt.
 *
 * @param c the #cell.
 * @param s the #scheduler.
 *
 * @return 1 If the space needs rebuilding. 0 otherwise.
 */
int cell_unskip_gravity_tasks(struct cell *c, struct scheduler *s) {

  struct engine *e = s->space->e;
  const int nodeID = e->nodeID;
  int rebuild = 0;

  /* Un-skip the gravity tasks involved with this cell. */
  for (struct link *l = c->grav.grav; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const int ci_active = cell_is_active_gravity(ci, e);
    const int cj_active = (cj != NULL) ? cell_is_active_gravity(cj, e) : 0;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active && ci_nodeID == nodeID) ||
        (cj_active && cj_nodeID == nodeID)) {

      scheduler_activate(s, t);

      /* Set the drifting flags */
      if (t->type == task_type_self &&
          t->subtype == task_subtype_external_grav) {
        cell_activate_subcell_external_grav_tasks(ci, s);
      } else if (t->type == task_type_self && t->subtype == task_subtype_grav) {
        cell_activate_subcell_grav_tasks(ci, NULL, s);
      } else if (t->type == task_type_pair) {
        cell_activate_subcell_grav_tasks(ci, cj, s);
      } else if (t->type == task_type_grav_mm) {
#ifdef SWIFT_DEBUG_CHECKS
        error("Incorrectly linked M-M task!");
#endif
      }
    }

    if (t->type == task_type_pair) {

#ifdef WITH_MPI
      /* Activate the send/recv tasks. */
      if (ci_nodeID != nodeID) {

        /* If the local cell is active, receive data from the foreign cell. */
        if (cj_active) scheduler_activate(s, ci->mpi.grav.recv);

        /* If the foreign cell is active, we want its ti_end values. */
        if (ci_active) scheduler_activate(s, ci->mpi.recv_ti);

        /* Is the foreign cell active and will need stuff from us? */
        if (ci_active) {

          scheduler_activate_send(s, cj->mpi.grav.send, ci_nodeID);

          /* Drift the cell which will be sent at the level at which it is
             sent, i.e. drift the cell specified in the send task (l->t)
             itself. */
          cell_activate_drift_gpart(cj, s);
        }

        /* If the local cell is active, send its ti_end values. */
        if (cj_active) scheduler_activate_send(s, cj->mpi.send_ti, ci_nodeID);

      } else if (cj_nodeID != nodeID) {

        /* If the local cell is active, receive data from the foreign cell. */
        if (ci_active) scheduler_activate(s, cj->mpi.grav.recv);

        /* If the foreign cell is active, we want its ti_end values. */
        if (cj_active) scheduler_activate(s, cj->mpi.recv_ti);

        /* Is the foreign cell active and will need stuff from us? */
        if (cj_active) {

          scheduler_activate_send(s, ci->mpi.grav.send, cj_nodeID);

          /* Drift the cell which will be sent at the level at which it is
             sent, i.e. drift the cell specified in the send task (l->t)
             itself. */
          cell_activate_drift_gpart(ci, s);
        }

        /* If the local cell is active, send its ti_end values. */
        if (ci_active) scheduler_activate_send(s, ci->mpi.send_ti, cj_nodeID);
      }
#endif
    }
  }

  for (struct link *l = c->grav.mm; l != NULL; l = l->next) {

    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const int ci_active = cell_is_active_gravity_mm(ci, e);
    const int cj_active = cell_is_active_gravity_mm(cj, e);
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

#ifdef SWIFT_DEBUG_CHECKS
    if (t->type != task_type_grav_mm) error("Incorrectly linked gravity task!");
#endif

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active && ci_nodeID == nodeID) ||
        (cj_active && cj_nodeID == nodeID)) {

      scheduler_activate(s, t);
    }
  }

  /* Unskip all the other task types. */
  if (c->nodeID == nodeID && cell_is_active_gravity(c, e)) {

    if (c->grav.init != NULL) scheduler_activate(s, c->grav.init);
    if (c->grav.init_out != NULL) scheduler_activate(s, c->grav.init_out);
    if (c->kick1 != NULL) scheduler_activate(s, c->kick1);
    if (c->kick2 != NULL) scheduler_activate(s, c->kick2);
    if (c->timestep != NULL) scheduler_activate(s, c->timestep);
    if (c->grav.down != NULL) scheduler_activate(s, c->grav.down);
    if (c->grav.down_in != NULL) scheduler_activate(s, c->grav.down_in);
    if (c->grav.mesh != NULL) scheduler_activate(s, c->grav.mesh);
    if (c->grav.long_range != NULL) scheduler_activate(s, c->grav.long_range);
    if (c->grav.end_force != NULL) scheduler_activate(s, c->grav.end_force);
    if (c->logger != NULL) scheduler_activate(s, c->logger);

    /* Subgrid tasks */
    if ((e->policy & engine_policy_cooling) && c->hydro.cooling != NULL)
      scheduler_activate(s, c->hydro.cooling);
    if ((e->policy & engine_policy_star_formation) &&
        c->hydro.star_formation != NULL)
      scheduler_activate(s, c->hydro.star_formation);
  }

  return rebuild;
}

/**
 * @brief Un-skips all the stars tasks associated with a given cell and checks
 * if the space needs to be rebuilt.
 *
 * @param c the #cell.
 * @param s the #scheduler.
 *
 * @return 1 If the space needs rebuilding. 0 otherwise.
 */
int cell_unskip_stars_tasks(struct cell *c, struct scheduler *s) {

  struct engine *e = s->space->e;
  const int with_feedback = (e->policy & engine_policy_feedback);
  const int nodeID = e->nodeID;
  int rebuild = 0;

  if (!with_feedback && c->stars.drift != NULL && cell_is_active_stars(c, e)) {
    cell_activate_drift_spart(c, s);
  }

  /* Un-skip the density tasks involved with this cell. */
  for (struct link *l = c->stars.density; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const int ci_active = cell_is_active_stars(ci, e);
    const int cj_active = (cj != NULL) ? cell_is_active_stars(cj, e) : 0;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    /* Activate the drifts */
    if (t->type == task_type_self && ci_active) {
      cell_activate_drift_part(ci, s);
      cell_activate_drift_spart(ci, s);
    }

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active || cj_active) &&
        (ci_nodeID == nodeID || cj_nodeID == nodeID)) {

      scheduler_activate(s, t);

      if (t->type == task_type_pair) {

        /* Do ci */
        if (ci_active) {
          /* stars for ci */
          atomic_or(&ci->stars.requires_sorts, 1 << t->flags);
          ci->stars.dx_max_sort_old = ci->stars.dx_max_sort;

          /* hydro for cj */
          atomic_or(&cj->hydro.requires_sorts, 1 << t->flags);
          cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;

          /* Activate the drift tasks. */
          if (ci_nodeID == nodeID) cell_activate_drift_spart(ci, s);
          if (cj_nodeID == nodeID) cell_activate_drift_part(cj, s);

          /* Check the sorts and activate them if needed. */
          cell_activate_stars_sorts(ci, t->flags, s);
          cell_activate_hydro_sorts(cj, t->flags, s);
        }

        /* Do cj */
        if (cj_active) {
          /* hydro for ci */
          atomic_or(&ci->hydro.requires_sorts, 1 << t->flags);
          ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;

          /* stars for cj */
          atomic_or(&cj->stars.requires_sorts, 1 << t->flags);
          cj->stars.dx_max_sort_old = cj->stars.dx_max_sort;

          /* Activate the drift tasks. */
          if (cj_nodeID == nodeID) cell_activate_drift_spart(cj, s);
          if (ci_nodeID == nodeID) cell_activate_drift_part(ci, s);

          /* Check the sorts and activate them if needed. */
          cell_activate_hydro_sorts(ci, t->flags, s);
          cell_activate_stars_sorts(cj, t->flags, s);
        }
      }

      else if (t->type == task_type_sub_pair || t->type == task_type_sub_self) {
        cell_activate_subcell_stars_tasks(ci, cj, s);
      }
    }

    /* Only interested in pair interactions as of here. */
    if (t->type == task_type_pair || t->type == task_type_sub_pair) {

      /* Check whether there was too much particle motion, i.e. the
         cell neighbour conditions were violated. */
      if (cell_need_rebuild_for_stars_pair(ci, cj)) rebuild = 1;
      if (cell_need_rebuild_for_stars_pair(cj, ci)) rebuild = 1;

#ifdef WITH_MPI
      /* Activate the send/recv tasks. */
      if (ci_nodeID != nodeID) {

        if (cj_active) {
          scheduler_activate(s, ci->mpi.hydro.recv_xv);
          scheduler_activate(s, ci->mpi.hydro.recv_rho);

          /* If the local cell is active, more stuff will be needed. */
          scheduler_activate_send(s, cj->mpi.stars.send, ci_nodeID);
          cell_activate_drift_spart(cj, s);

          /* If the local cell is active, send its ti_end values. */
          scheduler_activate_send(s, cj->mpi.send_ti, ci_nodeID);
        }

        if (ci_active) {
          scheduler_activate(s, ci->mpi.stars.recv);

          /* If the foreign cell is active, we want its ti_end values. */
          scheduler_activate(s, ci->mpi.recv_ti);

          /* Is the foreign cell active and will need stuff from us? */
          scheduler_activate_send(s, cj->mpi.hydro.send_xv, ci_nodeID);
          scheduler_activate_send(s, cj->mpi.hydro.send_rho, ci_nodeID);

          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(cj, s);
        }

      } else if (cj_nodeID != nodeID) {

        /* If the local cell is active, receive data from the foreign cell. */
        if (ci_active) {
          scheduler_activate(s, cj->mpi.hydro.recv_xv);
          scheduler_activate(s, cj->mpi.hydro.recv_rho);

          /* If the local cell is active, more stuff will be needed. */
          scheduler_activate_send(s, ci->mpi.stars.send, cj_nodeID);
          cell_activate_drift_spart(ci, s);

          /* If the local cell is active, send its ti_end values. */
          scheduler_activate_send(s, ci->mpi.send_ti, cj_nodeID);
        }

        if (cj_active) {
          scheduler_activate(s, cj->mpi.stars.recv);

          /* If the foreign cell is active, we want its ti_end values. */
          scheduler_activate(s, cj->mpi.recv_ti);

          /* Is the foreign cell active and will need stuff from us? */
          scheduler_activate_send(s, ci->mpi.hydro.send_xv, cj_nodeID);
          scheduler_activate_send(s, ci->mpi.hydro.send_rho, cj_nodeID);

          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(ci, s);
        }
      }
#endif
    }
  }

  /* Un-skip the feedback tasks involved with this cell. */
  for (struct link *l = c->stars.feedback; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const int ci_active = cell_is_active_stars(ci, e);
    const int cj_active = (cj != NULL) ? cell_is_active_stars(cj, e) : 0;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    if (t->type == task_type_self && ci_active) {
      scheduler_activate(s, t);
    }

    else if (t->type == task_type_sub_self && ci_active) {
      scheduler_activate(s, t);
    }

    else if (t->type == task_type_pair || t->type == task_type_sub_pair) {

      /* We only want to activate the task if the cell is active and is
         going to update some gas on the *local* node */
      if ((ci_nodeID == nodeID && cj_nodeID == nodeID) &&
          (ci_active || cj_active)) {

        scheduler_activate(s, t);

      } else if ((ci_nodeID == nodeID && cj_nodeID != nodeID) && (cj_active)) {

        scheduler_activate(s, t);

      } else if ((ci_nodeID != nodeID && cj_nodeID == nodeID) && (ci_active)) {

        scheduler_activate(s, t);
      }
    }

    /* Nothing more to do here, all drifts and sorts activated above */
  }

  /* Unskip all the other task types. */
  if (c->nodeID == nodeID && cell_is_active_stars(c, e)) {

    if (c->stars.ghost != NULL) scheduler_activate(s, c->stars.ghost);
    if (c->stars.stars_in != NULL) scheduler_activate(s, c->stars.stars_in);
    if (c->stars.stars_out != NULL) scheduler_activate(s, c->stars.stars_out);
    if (c->logger != NULL) scheduler_activate(s, c->logger);
  }

  return rebuild;
}

/**
 * @brief Set the super-cell pointers for all cells in a hierarchy.
 *
 * @param c The top-level #cell to play with.
 * @param super Pointer to the deepest cell with tasks in this part of the tree.
 * @param with_hydro Are we running with hydrodynamics on?
 * @param with_grav Are we running with gravity on?
 */
void cell_set_super(struct cell *c, struct cell *super, const int with_hydro,
                    const int with_grav) {

  /* Are we in a cell which is either the hydro or gravity super? */
  if (super == NULL && ((with_hydro && c->hydro.super != NULL) ||
                        (with_grav && c->grav.super != NULL)))
    super = c;

  /* Set the super-cell */
  c->super = super;

  /* Recurse */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        cell_set_super(c->progeny[k], super, with_hydro, with_grav);
}

/**
 * @brief Set the super-cell pointers for all cells in a hierarchy.
 *
 * @param c The top-level #cell to play with.
 * @param super_hydro Pointer to the deepest cell with tasks in this part of the
 * tree.
 */
void cell_set_super_hydro(struct cell *c, struct cell *super_hydro) {

  /* Are we in a cell with some kind of self/pair task ? */
  if (super_hydro == NULL && c->hydro.density != NULL) super_hydro = c;

  /* Set the super-cell */
  c->hydro.super = super_hydro;

  /* Recurse */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        cell_set_super_hydro(c->progeny[k], super_hydro);
}

/**
 * @brief Set the super-cell pointers for all cells in a hierarchy.
 *
 * @param c The top-level #cell to play with.
 * @param super_gravity Pointer to the deepest cell with tasks in this part of
 * the tree.
 */
void cell_set_super_gravity(struct cell *c, struct cell *super_gravity) {

  /* Are we in a cell with some kind of self/pair task ? */
  if (super_gravity == NULL && (c->grav.grav != NULL || c->grav.mm != NULL))
    super_gravity = c;

  /* Set the super-cell */
  c->grav.super = super_gravity;

  /* Recurse */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        cell_set_super_gravity(c->progeny[k], super_gravity);
}

/**
 * @brief Mapper function to set the super pointer of the cells.
 *
 * @param map_data The top-level cells.
 * @param num_elements The number of top-level cells.
 * @param extra_data Unused parameter.
 */
void cell_set_super_mapper(void *map_data, int num_elements, void *extra_data) {

  const struct engine *e = (const struct engine *)extra_data;

  const int with_hydro = (e->policy & engine_policy_hydro);
  const int with_grav = (e->policy & engine_policy_self_gravity) ||
                        (e->policy & engine_policy_external_gravity);

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &((struct cell *)map_data)[ind];

    /* All top-level cells get an MPI tag. */
#ifdef WITH_MPI
    cell_ensure_tagged(c);
#endif

    /* Super-pointer for hydro */
    if (with_hydro) cell_set_super_hydro(c, NULL);

    /* Super-pointer for gravity */
    if (with_grav) cell_set_super_gravity(c, NULL);

    /* Super-pointer for common operations */
    cell_set_super(c, NULL, with_hydro, with_grav);
  }
}

/**
 * @brief Does this cell or any of its children have any task ?
 *
 * We use the timestep-related tasks to probe this as these always
 * exist in a cell hierarchy that has any kind of task.
 *
 * @param c The #cell to probe.
 */
int cell_has_tasks(struct cell *c) {

#ifdef WITH_MPI
  if (c->timestep != NULL || c->mpi.recv_ti != NULL) return 1;
#else
  if (c->timestep != NULL) return 1;
#endif

  if (c->split) {
    int count = 0;
    for (int k = 0; k < 8; ++k)
      if (c->progeny[k] != NULL) count += cell_has_tasks(c->progeny[k]);
    return count;
  } else {
    return 0;
  }
}

/**
 * @brief Recursively drifts the #part in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 * @param force Drift the particles irrespective of the #cell flags.
 */
void cell_drift_part(struct cell *c, const struct engine *e, int force) {

  const int periodic = e->s->periodic;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const float hydro_h_max = e->hydro_properties->h_max;
  const float hydro_h_min = e->hydro_properties->h_min;
  const integertime_t ti_old_part = c->hydro.ti_old_part;
  const integertime_t ti_current = e->ti_current;
  struct part *const parts = c->hydro.parts;
  struct xpart *const xparts = c->hydro.xparts;

  float dx_max = 0.f, dx2_max = 0.f;
  float dx_max_sort = 0.0f, dx2_max_sort = 0.f;
  float cell_h_max = 0.f;

  /* Drift irrespective of cell flags? */
  force |= c->hydro.do_drift;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_part) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->hydro.count == 0) {

    /* Clear the drift flags. */
    c->hydro.do_drift = 0;
    c->hydro.do_sub_drift = 0;

    /* Update the time of the last drift */
    c->hydro.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || c->hydro.do_sub_drift)) {

    /* Loop over the progeny and collect their data. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];

        /* Collect */
        cell_drift_part(cp, e, force);

        /* Update */
        dx_max = max(dx_max, cp->hydro.dx_max_part);
        dx_max_sort = max(dx_max_sort, cp->hydro.dx_max_sort);
        cell_h_max = max(cell_h_max, cp->hydro.h_max);
      }
    }

    /* Store the values */
    c->hydro.h_max = cell_h_max;
    c->hydro.dx_max_part = dx_max;
    c->hydro.dx_max_sort = dx_max_sort;

    /* Update the time of the last drift */
    c->hydro.ti_old_part = ti_current;

  } else if (!c->split && force && ti_current > ti_old_part) {

    /* Drift from the last time the cell was drifted to the current time */
    double dt_drift, dt_kick_grav, dt_kick_hydro, dt_therm;
    if (with_cosmology) {
      dt_drift =
          cosmology_get_drift_factor(e->cosmology, ti_old_part, ti_current);
      dt_kick_grav =
          cosmology_get_grav_kick_factor(e->cosmology, ti_old_part, ti_current);
      dt_kick_hydro = cosmology_get_hydro_kick_factor(e->cosmology, ti_old_part,
                                                      ti_current);
      dt_therm = cosmology_get_therm_kick_factor(e->cosmology, ti_old_part,
                                                 ti_current);
    } else {
      dt_drift = (ti_current - ti_old_part) * e->time_base;
      dt_kick_grav = (ti_current - ti_old_part) * e->time_base;
      dt_kick_hydro = (ti_current - ti_old_part) * e->time_base;
      dt_therm = (ti_current - ti_old_part) * e->time_base;
    }

    /* Loop over all the gas particles in the cell */
    const size_t nr_parts = c->hydro.count;
    for (size_t k = 0; k < nr_parts; k++) {

      /* Get a handle on the part. */
      struct part *const p = &parts[k];
      struct xpart *const xp = &xparts[k];

      /* Ignore inhibited particles */
      if (part_is_inhibited(p, e)) continue;

      /* Drift... */
      drift_part(p, xp, dt_drift, dt_kick_hydro, dt_kick_grav, dt_therm,
                 ti_old_part, ti_current);

      /* Update the tracers properties */
      tracers_after_drift(p, xp, e->internal_units, e->physical_constants,
                          with_cosmology, e->cosmology, e->hydro_properties,
                          e->cooling_func, e->time);

#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure the particle does not drift by more than a box length. */
      if (fabs(xp->v_full[0] * dt_drift) > e->s->dim[0] ||
          fabs(xp->v_full[1] * dt_drift) > e->s->dim[1] ||
          fabs(xp->v_full[2] * dt_drift) > e->s->dim[2]) {
        error("Particle drifts by more than a box length!");
      }
#endif

      /* In non-periodic BC runs, remove particles that crossed the border */
      if (!periodic) {

        /* Did the particle leave the box?  */
        if ((p->x[0] > dim[0]) || (p->x[0] < 0.) ||  // x
            (p->x[1] > dim[1]) || (p->x[1] < 0.) ||  // y
            (p->x[2] > dim[2]) || (p->x[2] < 0.)) {  // z

          /* One last action before death? */
          hydro_remove_part(p, xp);

          /* Remove the particle entirely */
          struct gpart *gp = p->gpart;
          cell_remove_part(e, c, p, xp);

          /* and it's gravity friend */
          if (gp != NULL) cell_remove_gpart(e, c, gp);

          continue;
        }
      }

      /* Limit h to within the allowed range */
      p->h = min(p->h, hydro_h_max);
      p->h = max(p->h, hydro_h_min);

      /* Compute (square of) motion since last cell construction */
      const float dx2 = xp->x_diff[0] * xp->x_diff[0] +
                        xp->x_diff[1] * xp->x_diff[1] +
                        xp->x_diff[2] * xp->x_diff[2];
      dx2_max = max(dx2_max, dx2);
      const float dx2_sort = xp->x_diff_sort[0] * xp->x_diff_sort[0] +
                             xp->x_diff_sort[1] * xp->x_diff_sort[1] +
                             xp->x_diff_sort[2] * xp->x_diff_sort[2];
      dx2_max_sort = max(dx2_max_sort, dx2_sort);

      /* Maximal smoothing length */
      cell_h_max = max(cell_h_max, p->h);

      /* Get ready for a density calculation */
      if (part_is_active(p, e)) {
        hydro_init_part(p, &e->s->hs);
        chemistry_init_part(p, e->chemistry);
        tracers_after_init(p, xp, e->internal_units, e->physical_constants,
                           with_cosmology, e->cosmology, e->hydro_properties,
                           e->cooling_func, e->time);
      }
    }

    /* Now, get the maximal particle motion from its square */
    dx_max = sqrtf(dx2_max);
    dx_max_sort = sqrtf(dx2_max_sort);

    /* Store the values */
    c->hydro.h_max = cell_h_max;
    c->hydro.dx_max_part = dx_max;
    c->hydro.dx_max_sort = dx_max_sort;

    /* Update the time of the last drift */
    c->hydro.ti_old_part = ti_current;
  }

  /* Clear the drift flags. */
  c->hydro.do_drift = 0;
  c->hydro.do_sub_drift = 0;
}

/**
 * @brief Recursively drifts the #gpart in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 * @param force Drift the particles irrespective of the #cell flags.
 */
void cell_drift_gpart(struct cell *c, const struct engine *e, int force) {

  const int periodic = e->s->periodic;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const integertime_t ti_old_gpart = c->grav.ti_old_part;
  const integertime_t ti_current = e->ti_current;
  struct gpart *const gparts = c->grav.parts;

  /* Drift irrespective of cell flags? */
  force |= c->grav.do_drift;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_gpart) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->grav.count == 0) {

    /* Clear the drift flags. */
    c->grav.do_drift = 0;
    c->grav.do_sub_drift = 0;

    /* Update the time of the last drift */
    c->grav.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || c->grav.do_sub_drift)) {

    /* Loop over the progeny and collect their data. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];

        /* Recurse */
        cell_drift_gpart(cp, e, force);
      }
    }

    /* Update the time of the last drift */
    c->grav.ti_old_part = ti_current;

  } else if (!c->split && force && ti_current > ti_old_gpart) {

    /* Drift from the last time the cell was drifted to the current time */
    double dt_drift;
    if (with_cosmology) {
      dt_drift =
          cosmology_get_drift_factor(e->cosmology, ti_old_gpart, ti_current);
    } else {
      dt_drift = (ti_current - ti_old_gpart) * e->time_base;
    }

    /* Loop over all the g-particles in the cell */
    const size_t nr_gparts = c->grav.count;
    for (size_t k = 0; k < nr_gparts; k++) {

      /* Get a handle on the gpart. */
      struct gpart *const gp = &gparts[k];

      /* Ignore inhibited particles */
      if (gpart_is_inhibited(gp, e)) continue;

      /* Drift... */
      drift_gpart(gp, dt_drift, ti_old_gpart, ti_current);

#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure the particle does not drift by more than a box length. */
      if (fabs(gp->v_full[0] * dt_drift) > e->s->dim[0] ||
          fabs(gp->v_full[1] * dt_drift) > e->s->dim[1] ||
          fabs(gp->v_full[2] * dt_drift) > e->s->dim[2]) {
        error("Particle drifts by more than a box length!");
      }
#endif

      /* In non-periodic BC runs, remove particles that crossed the border */
      if (!periodic) {

        /* Did the particle leave the box?  */
        if ((gp->x[0] > dim[0]) || (gp->x[0] < 0.) ||  // x
            (gp->x[1] > dim[1]) || (gp->x[1] < 0.) ||  // y
            (gp->x[2] > dim[2]) || (gp->x[2] < 0.)) {  // z

          /* Remove the particle entirely */
          if (gp->type == swift_type_dark_matter) cell_remove_gpart(e, c, gp);

          continue;
        }
      }

      /* Init gravity force fields. */
      if (gpart_is_active(gp, e)) {
        gravity_init_gpart(gp);
      }
    }

    /* Update the time of the last drift */
    c->grav.ti_old_part = ti_current;
  }

  /* Clear the drift flags. */
  c->grav.do_drift = 0;
  c->grav.do_sub_drift = 0;
}

/**
 * @brief Recursively drifts the #spart in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 * @param force Drift the particles irrespective of the #cell flags.
 */
void cell_drift_spart(struct cell *c, const struct engine *e, int force) {

  const int periodic = e->s->periodic;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const float stars_h_max = e->hydro_properties->h_max;
  const float stars_h_min = e->hydro_properties->h_min;
  const integertime_t ti_old_spart = c->stars.ti_old_part;
  const integertime_t ti_current = e->ti_current;
  struct spart *const sparts = c->stars.parts;

  float dx_max = 0.f, dx2_max = 0.f;
  float dx_max_sort = 0.0f, dx2_max_sort = 0.f;
  float cell_h_max = 0.f;

  /* Drift irrespective of cell flags? */
  force |= c->stars.do_drift;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_spart) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->stars.count == 0) {

    /* Clear the drift flags. */
    c->stars.do_drift = 0;
    c->stars.do_sub_drift = 0;

    /* Update the time of the last drift */
    c->stars.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || c->stars.do_sub_drift)) {

    /* Loop over the progeny and collect their data. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];

        /* Recurse */
        cell_drift_spart(cp, e, force);

        /* Update */
        dx_max = max(dx_max, cp->stars.dx_max_part);
        dx_max_sort = max(dx_max_sort, cp->stars.dx_max_sort);
        cell_h_max = max(cell_h_max, cp->stars.h_max);
      }
    }

    /* Store the values */
    c->stars.h_max = cell_h_max;
    c->stars.dx_max_part = dx_max;
    c->stars.dx_max_sort = dx_max_sort;

    /* Update the time of the last drift */
    c->stars.ti_old_part = ti_current;

  } else if (!c->split && force && ti_current > ti_old_spart) {

    /* Drift from the last time the cell was drifted to the current time */
    double dt_drift;
    if (with_cosmology) {
      dt_drift =
          cosmology_get_drift_factor(e->cosmology, ti_old_spart, ti_current);
    } else {
      dt_drift = (ti_current - ti_old_spart) * e->time_base;
    }

    /* Loop over all the star particles in the cell */
    const size_t nr_sparts = c->stars.count;
    for (size_t k = 0; k < nr_sparts; k++) {

      /* Get a handle on the spart. */
      struct spart *const sp = &sparts[k];

      /* Ignore inhibited particles */
      if (spart_is_inhibited(sp, e)) continue;

      /* Drift... */
      drift_spart(sp, dt_drift, ti_old_spart, ti_current);

#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure the particle does not drift by more than a box length. */
      if (fabs(sp->v[0] * dt_drift) > e->s->dim[0] ||
          fabs(sp->v[1] * dt_drift) > e->s->dim[1] ||
          fabs(sp->v[2] * dt_drift) > e->s->dim[2]) {
        error("Particle drifts by more than a box length!");
      }
#endif

      /* In non-periodic BC runs, remove particles that crossed the border */
      if (!periodic) {

        /* Did the particle leave the box?  */
        if ((sp->x[0] > dim[0]) || (sp->x[0] < 0.) ||  // x
            (sp->x[1] > dim[1]) || (sp->x[1] < 0.) ||  // y
            (sp->x[2] > dim[2]) || (sp->x[2] < 0.)) {  // z

          /* Remove the particle entirely */
          struct gpart *gp = sp->gpart;
          cell_remove_spart(e, c, sp);

          /* and it's gravity friend */
          cell_remove_gpart(e, c, gp);

          continue;
        }
      }

      /* Limit h to within the allowed range */
      sp->h = min(sp->h, stars_h_max);
      sp->h = max(sp->h, stars_h_min);

      /* Compute (square of) motion since last cell construction */
      const float dx2 = sp->x_diff[0] * sp->x_diff[0] +
                        sp->x_diff[1] * sp->x_diff[1] +
                        sp->x_diff[2] * sp->x_diff[2];
      dx2_max = max(dx2_max, dx2);

      const float dx2_sort = sp->x_diff_sort[0] * sp->x_diff_sort[0] +
                             sp->x_diff_sort[1] * sp->x_diff_sort[1] +
                             sp->x_diff_sort[2] * sp->x_diff_sort[2];

      dx2_max_sort = max(dx2_max_sort, dx2_sort);

      /* Maximal smoothing length */
      cell_h_max = max(cell_h_max, sp->h);

      /* Get ready for a density calculation */
      if (spart_is_active(sp, e)) {
        stars_init_spart(sp);
      }
    }

    /* Now, get the maximal particle motion from its square */
    dx_max = sqrtf(dx2_max);
    dx_max_sort = sqrtf(dx2_max_sort);

    /* Store the values */
    c->stars.h_max = cell_h_max;
    c->stars.dx_max_part = dx_max;
    c->stars.dx_max_sort = dx_max_sort;

    /* Update the time of the last drift */
    c->stars.ti_old_part = ti_current;
  }

  /* Clear the drift flags. */
  c->stars.do_drift = 0;
  c->stars.do_sub_drift = 0;
}

/**
 * @brief Recursively drifts all multipoles in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 */
void cell_drift_all_multipoles(struct cell *c, const struct engine *e) {

  const integertime_t ti_old_multipole = c->grav.ti_old_multipole;
  const integertime_t ti_current = e->ti_current;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_multipole) error("Attempt to drift to the past");
#endif

  /* Drift from the last time the cell was drifted to the current time */
  double dt_drift;
  if (e->policy & engine_policy_cosmology)
    dt_drift =
        cosmology_get_drift_factor(e->cosmology, ti_old_multipole, ti_current);
  else
    dt_drift = (ti_current - ti_old_multipole) * e->time_base;

  /* Drift the multipole */
  if (ti_current > ti_old_multipole) gravity_drift(c->grav.multipole, dt_drift);

  /* Are we not in a leaf ? */
  if (c->split) {

    /* Loop over the progeny and recurse. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_drift_all_multipoles(c->progeny[k], e);
  }

  /* Update the time of the last drift */
  c->grav.ti_old_multipole = ti_current;
}

/**
 * @brief Drifts the multipole of a cell to the current time.
 *
 * Only drifts the multipole at this level. Multipoles deeper in the
 * tree are not updated.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 */
void cell_drift_multipole(struct cell *c, const struct engine *e) {

  const integertime_t ti_old_multipole = c->grav.ti_old_multipole;
  const integertime_t ti_current = e->ti_current;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_multipole) error("Attempt to drift to the past");
#endif

  /* Drift from the last time the cell was drifted to the current time */
  double dt_drift;
  if (e->policy & engine_policy_cosmology)
    dt_drift =
        cosmology_get_drift_factor(e->cosmology, ti_old_multipole, ti_current);
  else
    dt_drift = (ti_current - ti_old_multipole) * e->time_base;

  if (ti_current > ti_old_multipole) gravity_drift(c->grav.multipole, dt_drift);

  /* Update the time of the last drift */
  c->grav.ti_old_multipole = ti_current;
}

/**
 * @brief Resets all the sorting properties for the stars in a given cell
 * hierarchy.
 *
 * @param c The #cell to clean.
 * @param is_super Is this a super-cell?
 */
void cell_clear_stars_sort_flags(struct cell *c, const int is_super) {

  /* Recurse if possible */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        cell_clear_stars_sort_flags(c->progeny[k], /*is_super=*/0);
  }

  /* Free the sorted array at the level where it was allocated */
  if (is_super) {

#ifdef SWIFT_DEBUG_CHECKS
    if (c != c->hydro.super) error("Cell is not a super-cell!!!");
#endif
    cell_free_stars_sorts(c);
  }

  /* Indicate that the cell is not sorted and cancel the pointer sorting arrays.
   */
  c->stars.sorted = 0;
  for (int i = 0; i < 13; i++) {
    c->stars.sort[i] = NULL;
  }
}

/**
 * @brief Recursively checks that all particles in a cell have a time-step
 */
void cell_check_timesteps(struct cell *c) {
#ifdef SWIFT_DEBUG_CHECKS

  if (c->hydro.ti_end_min == 0 && c->grav.ti_end_min == 0 &&
      c->stars.ti_end_min == 0 && c->nr_tasks > 0)
    error("Cell without assigned time-step");

  if (c->split) {
    for (int k = 0; k < 8; ++k)
      if (c->progeny[k] != NULL) cell_check_timesteps(c->progeny[k]);
  } else {

    if (c->nodeID == engine_rank)
      for (int i = 0; i < c->hydro.count; ++i)
        if (c->hydro.parts[i].time_bin == 0)
          error("Particle without assigned time-bin");
  }
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

void cell_check_spart_pos(const struct cell *c,
                          const struct spart *global_sparts) {

#ifdef SWIFT_DEBUG_CHECKS

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; ++k)
      if (c->progeny[k] != NULL)
        cell_check_spart_pos(c->progeny[k], global_sparts);
  }

  struct spart *sparts = c->stars.parts;
  const int count = c->stars.count;
  for (int i = 0; i < count; ++i) {

    const struct spart *sp = &sparts[i];
    if ((sp->x[0] < c->loc[0] / space_stretch) ||
        (sp->x[1] < c->loc[1] / space_stretch) ||
        (sp->x[2] < c->loc[2] / space_stretch) ||
        (sp->x[0] >= (c->loc[0] + c->width[0]) * space_stretch) ||
        (sp->x[1] >= (c->loc[1] + c->width[1]) * space_stretch) ||
        (sp->x[2] >= (c->loc[2] + c->width[2]) * space_stretch))
      error("spart not in its cell!");

    if (sp->time_bin != time_bin_not_created &&
        sp->time_bin != time_bin_inhibited) {

      const struct gpart *gp = sp->gpart;
      if (gp == NULL && sp->time_bin != time_bin_not_created)
        error("Unlinked spart!");

      if (&global_sparts[-gp->id_or_neg_offset] != sp)
        error("Incorrectly linked spart!");
    }
  }

#else
  error("Calling a degugging function outside debugging mode.");
#endif
}

/**
 * @brief Recursively update the pointer and counter for #spart after the
 * addition of a new particle.
 *
 * @param c The cell we are working on.
 * @param progeny_list The list of the progeny index at each level for the
 * leaf-cell where the particle was added.
 * @param main_branch Are we in a cell directly above the leaf where the new
 * particle was added?
 */
void cell_recursively_shift_sparts(struct cell *c,
                                   const int progeny_list[space_cell_maxdepth],
                                   const int main_branch) {
  if (c->split) {

    /* No need to recurse in progenies located before the insestion point */
    const int first_progeny = main_branch ? progeny_list[(int)c->depth] : 0;

    for (int k = first_progeny; k < 8; ++k) {

      if (c->progeny[k] != NULL)
        cell_recursively_shift_sparts(c->progeny[k], progeny_list,
                                      main_branch && (k == first_progeny));
    }
  }

  /* When directly above the leaf with the new particle: increase the particle
   * count */
  /* When after the leaf with the new particle: shift by one position */
  if (main_branch)
    c->stars.count++;
  else
    c->stars.parts++;
}

/**
 * @brief "Add" a #spart in a given #cell.
 *
 * This function will a a #spart at the start of the current cell's array by
 * shifting all the #spart in the top-level cell by one position. All the
 * pointers and cell counts are updated accordingly.
 *
 * @param e The #engine.
 * @param c The leaf-cell in which to add the #spart.
 *
 * @return A pointer to the newly added #spart. The spart has a been zeroed and
 * given a position within the cell as well as set to the minimal active time
 * bin.
 */
struct spart *cell_add_spart(struct engine *e, struct cell *const c) {

  /* Perform some basic consitency checks */
  if (c->nodeID != engine_rank) error("Adding spart on a foreign node");
  if (c->grav.ti_old_part != e->ti_current) error("Undrifted cell!");
  if (c->split) error("Addition of spart performed above the leaf level");

  /* Progeny number at each level */
  int progeny[space_cell_maxdepth];
#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < space_cell_maxdepth; ++i) progeny[i] = -1;
#endif

  /* Get the top-level this leaf cell is in and compute the progeny indices at
     each level */
  struct cell *top = c;
  while (top->parent != NULL) {
    for (int k = 0; k < 8; ++k) {
      if (top->parent->progeny[k] == top) {
        progeny[(int)top->parent->depth] = k;
      }
    }
    top = top->parent;
  }

  /* Lock the top-level cell as we are going to operate on it */
  lock_lock(&top->stars.star_formation_lock);

  /* Are there any extra particles left? */
  if (top->stars.count == top->stars.count_total - 1) {
    message("We ran out of star particles!");
    atomic_inc(&e->forcerebuild);
    return NULL;
  }

  /* Number of particles to shift in order to get a free space. */
  const size_t n_copy = &top->stars.parts[top->stars.count] - c->stars.parts;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->stars.parts + n_copy > top->stars.parts + top->stars.count)
    error("Copying beyond the allowed range");
#endif

  if (n_copy > 0) {

    // MATTHIEU: This can be improved. We don't need to copy everything, just
    // need to swap a few particles.
    memmove(&c->stars.parts[1], &c->stars.parts[0],
            n_copy * sizeof(struct spart));

    /* Update the gpart->spart links (shift by 1) */
    for (size_t i = 0; i < n_copy; ++i) {
#ifdef SWIFT_DEBUG_CHECKS
      if (c->stars.parts[i + 1].gpart == NULL) {
        error("Incorrectly linked spart!");
      }
#endif
      c->stars.parts[i + 1].gpart->id_or_neg_offset--;
    }
  }

  /* Recursively shift all the stars to get a free spot at the start of the
   * current cell*/
  cell_recursively_shift_sparts(top, progeny, /* main_branch=*/1);

  /* Make sure the gravity will be recomputed for this particle in the next step
   */
  struct cell *top2 = c;
  while (top2->parent != NULL) {
    top2->grav.ti_end_min = e->ti_current;
    top2 = top2->parent;
  }
  top2->grav.ti_end_min = e->ti_current;

  /* Release the lock */
  if (lock_unlock(&top->stars.star_formation_lock) != 0)
    error("Failed to unlock the top-level cell.");

  /* We now have an empty spart as the first particle in that cell */
  struct spart *sp = &c->stars.parts[0];
  bzero(sp, sizeof(struct spart));

  /* Give it a decent position */
  sp->x[0] = c->loc[0] + 0.5 * c->width[0];
  sp->x[1] = c->loc[1] + 0.5 * c->width[1];
  sp->x[2] = c->loc[2] + 0.5 * c->width[2];

  /* Set it to the current time-bin */
  sp->time_bin = e->min_active_bin;

#ifdef SWIFT_DEBUG_CHECKS
  /* Specify it was drifted to this point */
  sp->ti_drift = e->ti_current;
#endif

  /* Register that we used one of the free slots. */
  const size_t one = 1;
  atomic_sub(&e->s->nr_extra_sparts, one);

  return sp;
}

/**
 * @brief "Remove" a gas particle from the calculation.
 *
 * The particle is inhibited and will officially be removed at the next rebuild.
 *
 * @param e The #engine running on this node.
 * @param c The #cell from which to remove the particle.
 * @param p The #part to remove.
 * @param xp The extended data of the particle to remove.
 */
void cell_remove_part(const struct engine *e, struct cell *c, struct part *p,
                      struct xpart *xp) {

  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  /* Mark the particle as inhibited */
  p->time_bin = time_bin_inhibited;

  /* Mark the gpart as inhibited and stand-alone */
  if (p->gpart) {
    p->gpart->time_bin = time_bin_inhibited;
    p->gpart->id_or_neg_offset = p->id;
    p->gpart->type = swift_type_dark_matter;
  }

  /* Un-link the part */
  p->gpart = NULL;
}

/**
 * @brief "Remove" a gravity particle from the calculation.
 *
 * The particle is inhibited and will officially be removed at the next rebuild.
 *
 * @param e The #engine running on this node.
 * @param c The #cell from which to remove the particle.
 * @param gp The #gpart to remove.
 */
void cell_remove_gpart(const struct engine *e, struct cell *c,
                       struct gpart *gp) {

  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  if (gp->type != swift_type_dark_matter)
    error("Trying to remove a non-dark matter gpart.");

  /* Mark the particle as inhibited */
  gp->time_bin = time_bin_inhibited;
}

/**
 * @brief "Remove" a star particle from the calculation.
 *
 * The particle is inhibited and will officially be removed at the next rebuild.
 *
 * @param e The #engine running on this node.
 * @param c The #cell from which to remove the particle.
 * @param sp The #spart to remove.
 */
void cell_remove_spart(const struct engine *e, struct cell *c,
                       struct spart *sp) {

  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  /* Mark the particle as inhibited and stand-alone */
  sp->time_bin = time_bin_inhibited;
  if (sp->gpart) {
    sp->gpart->time_bin = time_bin_inhibited;
    sp->gpart->id_or_neg_offset = sp->id;
    sp->gpart->type = swift_type_dark_matter;
  }

  /* Un-link the spart */
  sp->gpart = NULL;
}

/**
 * @brief "Remove" a gas particle from the calculation and convert its gpart
 * friend to a dark matter particle.
 *
 * Note that the #part is not destroyed. The pointer is still valid
 * after this call and the properties of the #part are not altered
 * apart from the time-bin and #gpart pointer.
 * The particle is inhibited and will officially be removed at the next rebuild.
 *
 * @param e The #engine running on this node.
 * @param c The #cell from which to remove the particle.
 * @param p The #part to remove.
 * @param xp The extended data of the particle to remove.
 *
 * @return Pointer to the #gpart the #part has become. It carries the
 * ID of the #part and has a dark matter type.
 */
struct gpart *cell_convert_part_to_gpart(const struct engine *e, struct cell *c,
                                         struct part *p, struct xpart *xp) {

  /* Quick cross-checks */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  if (p->gpart == NULL)
    error("Trying to convert part without gpart friend to dark matter!");

  /* Get a handle */
  struct gpart *gp = p->gpart;

  /* Mark the particle as inhibited */
  p->time_bin = time_bin_inhibited;

  /* Un-link the part */
  p->gpart = NULL;

  /* Mark the gpart as dark matter */
  gp->type = swift_type_dark_matter;
  gp->id_or_neg_offset = p->id;

#ifdef SWIFT_DEBUG_CHECKS
  gp->ti_kick = p->ti_kick;
#endif

  return gp;
}

/**
 * @brief "Remove" a spart particle from the calculation and convert its gpart
 * friend to a dark matter particle.
 *
 * Note that the #spart is not destroyed. The pointer is still valid
 * after this call and the properties of the #spart are not altered
 * apart from the time-bin and #gpart pointer.
 * The particle is inhibited and will officially be removed at the next rebuild.
 *
 * @param e The #engine running on this node.
 * @param c The #cell from which to remove the particle.
 * @param sp The #spart to remove.
 *
 * @return Pointer to the #gpart the #spart has become. It carries the
 * ID of the #spart and has a dark matter type.
 */
struct gpart *cell_convert_spart_to_gpart(const struct engine *e,
                                          struct cell *c, struct spart *sp) {

  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  if (sp->gpart == NULL)
    error("Trying to convert spart without gpart friend to dark matter!");

  /* Get a handle */
  struct gpart *gp = sp->gpart;

  /* Mark the particle as inhibited */
  sp->time_bin = time_bin_inhibited;

  /* Un-link the spart */
  sp->gpart = NULL;

  /* Mark the gpart as dark matter */
  gp->type = swift_type_dark_matter;
  gp->id_or_neg_offset = sp->id;

#ifdef SWIFT_DEBUG_CHECKS
  gp->ti_kick = sp->ti_kick;
#endif

  return gp;
}

/**
 * @brief "Remove" a #part from a #cell and replace it with a #spart
 * connected to the same #gpart.
 *
 * Note that the #part is not destroyed. The pointer is still valid
 * after this call and the properties of the #part are not altered
 * apart from the time-bin and #gpart pointer.
 * The particle is inhibited and will officially be removed at the next rebuild.
 *
 * @param e The #engine.
 * @param c The #cell from which to remove the #part.
 * @param p The #part to remove (must be inside c).
 * @param xp The extended data of the #part.
 *
 * @return A fresh #spart with the same ID, position, velocity and
 * time-bin as the original #part.
 */
struct spart *cell_convert_part_to_spart(struct engine *e, struct cell *c,
                                         struct part *p, struct xpart *xp) {

  /* Quick cross-check */
  if (c->nodeID != e->nodeID)
    error("Can't remove a particle in a foreign cell.");

  if (p->gpart == NULL)
    error("Trying to convert part without gpart friend to star!");

  /* Create a fresh (empty) spart */
  struct spart *sp = cell_add_spart(e, c);

  /* Did we run out of free spart slots? */
  if (sp == NULL) return NULL;

  /* Destroy the gas particle and get it's gpart friend */
  struct gpart *gp = cell_convert_part_to_gpart(e, c, p, xp);

  /* Assign the ID back */
  sp->id = gp->id_or_neg_offset;
  gp->type = swift_type_stars;

  /* Re-link things */
  sp->gpart = gp;
  gp->id_or_neg_offset = -(sp - e->s->sparts);

  /* Synchronize clocks */
  gp->time_bin = sp->time_bin;

  /* Synchronize masses, positions and velocities */
  sp->mass = gp->mass;
  sp->x[0] = gp->x[0];
  sp->x[1] = gp->x[1];
  sp->x[2] = gp->x[2];
  sp->v[0] = gp->v_full[0];
  sp->v[1] = gp->v_full[1];
  sp->v[2] = gp->v_full[2];

#ifdef SWIFT_DEBUG_CHECKS
  sp->ti_kick = gp->ti_kick;
  gp->ti_drift = sp->ti_drift;
#endif

  /* Set a smoothing length */
  sp->h = max(c->stars.h_max, c->hydro.h_max);

  /* Here comes the Sun! */
  return sp;
}

/**
 * @brief Re-arrange the #part in a top-level cell such that all the extra ones
 * for on-the-fly creation are located at the end of the array.
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
}

/**
 * @brief Re-arrange the #spart in a top-level cell such that all the extra ones
 * for on-the-fly creation are located at the end of the array.
 *
 * @param c The #cell to sort.
 * @param sparts_offset The offset between the first #spart in the array and the
 * first #spart in the global array in the space structure (for re-linking).
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
}

/**
 * @brief Re-arrange the #gpart in a top-level cell such that all the extra ones
 * for on-the-fly creation are located at the end of the array.
 *
 * @param c The #cell to sort.
 * @param parts The global array of #part (for re-linking).
 * @param sparts The global array of #spart (for re-linking).
 */
void cell_reorder_extra_gparts(struct cell *c, struct part *parts,
                               struct spart *sparts) {

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
      memswap(&gparts[i], &gparts[first_not_extra], sizeof(struct gpart));
      if (gparts[i].type == swift_type_gas) {
        parts[-gparts[i].id_or_neg_offset].gpart = &gparts[i];
      } else if (gparts[i].type == swift_type_stars) {
        sparts[-gparts[i].id_or_neg_offset].gpart = &gparts[i];
      }
    }
  }
}

/**
 * @brief Can we use the MM interactions fo a given pair of cells?
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param e The #engine.
 * @param s The #space.
 */
int cell_can_use_pair_mm(const struct cell *ci, const struct cell *cj,
                         const struct engine *e, const struct space *s) {

  const double theta_crit2 = e->gravity_properties->theta_crit2;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  /* Recover the multipole information */
  const struct gravity_tensors *const multi_i = ci->grav.multipole;
  const struct gravity_tensors *const multi_j = cj->grav.multipole;

  /* Get the distance between the CoMs */
  double dx = multi_i->CoM[0] - multi_j->CoM[0];
  double dy = multi_i->CoM[1] - multi_j->CoM[1];
  double dz = multi_i->CoM[2] - multi_j->CoM[2];

  /* Apply BC */
  if (periodic) {
    dx = nearest(dx, dim[0]);
    dy = nearest(dy, dim[1]);
    dz = nearest(dz, dim[2]);
  }
  const double r2 = dx * dx + dy * dy + dz * dz;

  return gravity_M2L_accept(multi_i->r_max, multi_j->r_max, theta_crit2, r2);
}

/**
 * @brief Can we use the MM interactions fo a given pair of cells?
 *
 * This function uses the information gathered in the multipole at rebuild
 * time and not the current position and radius of the multipole.
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param e The #engine.
 * @param s The #space.
 */
int cell_can_use_pair_mm_rebuild(const struct cell *ci, const struct cell *cj,
                                 const struct engine *e,
                                 const struct space *s) {

  const double theta_crit2 = e->gravity_properties->theta_crit2;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  /* Recover the multipole information */
  const struct gravity_tensors *const multi_i = ci->grav.multipole;
  const struct gravity_tensors *const multi_j = cj->grav.multipole;

#ifdef SWIFT_DEBUG_CHECKS

  if (multi_i->CoM_rebuild[0] < ci->loc[0] ||
      multi_i->CoM_rebuild[0] > ci->loc[0] + ci->width[0])
    error("Invalid multipole position ci");
  if (multi_i->CoM_rebuild[1] < ci->loc[1] ||
      multi_i->CoM_rebuild[1] > ci->loc[1] + ci->width[1])
    error("Invalid multipole position ci");
  if (multi_i->CoM_rebuild[2] < ci->loc[2] ||
      multi_i->CoM_rebuild[2] > ci->loc[2] + ci->width[2])
    error("Invalid multipole position ci");

  if (multi_j->CoM_rebuild[0] < cj->loc[0] ||
      multi_j->CoM_rebuild[0] > cj->loc[0] + cj->width[0])
    error("Invalid multipole position cj");
  if (multi_j->CoM_rebuild[1] < cj->loc[1] ||
      multi_j->CoM_rebuild[1] > cj->loc[1] + cj->width[1])
    error("Invalid multipole position cj");
  if (multi_j->CoM_rebuild[2] < cj->loc[2] ||
      multi_j->CoM_rebuild[2] > cj->loc[2] + cj->width[2])
    error("Invalid multipole position cj");

#endif

  /* Get the distance between the CoMs */
  double dx = multi_i->CoM_rebuild[0] - multi_j->CoM_rebuild[0];
  double dy = multi_i->CoM_rebuild[1] - multi_j->CoM_rebuild[1];
  double dz = multi_i->CoM_rebuild[2] - multi_j->CoM_rebuild[2];

  /* Apply BC */
  if (periodic) {
    dx = nearest(dx, dim[0]);
    dy = nearest(dy, dim[1]);
    dz = nearest(dz, dim[2]);
  }
  const double r2 = dx * dx + dy * dy + dz * dz;

  return gravity_M2L_accept(multi_i->r_max_rebuild, multi_j->r_max_rebuild,
                            theta_crit2, r2);
}
