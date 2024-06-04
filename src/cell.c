/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include <config.h>

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
#include "engine.h"
#include "error.h"
#include "multipole.h"
#include "space.h"
#include "tools.h"

/* Global variables. */
int cell_next_tag = 0;

/** List of cell pairs for sub-cell recursion. For any sid, the entries in
 * this array contain the number of sub-cell pairs and the indices and sid
 * of the sub-cell pairs themselves. */
struct cell_split_pair cell_split_pairs[13] = {
    {1, /* (  1 ,  1 ,  1 ) */
     {{7, 0, 0}}},

    {4, /* (  1 ,  1 ,  0 ) */
     {{6, 0, 1}, {7, 1, 1}, {6, 1, 0}, {7, 0, 2}}},

    {1, /* (  1 ,  1 , -1 ) */
     {{6, 1, 2}}},

    {4, /* (  1 ,  0 ,  1 ) */
     {{5, 0, 3}, {7, 2, 3}, {5, 2, 0}, {7, 0, 6}}},

    {16, /* (  1 ,  0 ,  0 ) */
     {{4, 0, 4},
      {5, 0, 5},
      {6, 0, 7},
      {7, 0, 8},
      {4, 1, 3},
      {5, 1, 4},
      {6, 1, 6},
      {7, 1, 7},
      {4, 2, 1},
      {5, 2, 2},
      {6, 2, 4},
      {7, 2, 5},
      {4, 3, 0},
      {5, 3, 1},
      {6, 3, 3},
      {7, 3, 4}}},

    {4, /* (  1 ,  0 , -1 ) */
     {{4, 1, 5}, {6, 3, 5}, {4, 3, 2}, {6, 1, 8}}},

    {1, /* (  1 , -1 ,  1 ) */
     {{5, 2, 6}}},

    {4, /* (  1 , -1 ,  0 ) */
     {{4, 3, 6}, {5, 2, 8}, {4, 2, 7}, {5, 3, 7}}},

    {1, /* (  1 , -1 , -1 ) */
     {{4, 3, 8}}},

    {4, /* (  0 ,  1 ,  1 ) */
     {{3, 0, 9}, {7, 4, 9}, {3, 4, 0}, {7, 0, 8}}},

    {16, /* (  0 ,  1 ,  0 ) */
     {{2, 0, 10},
      {3, 0, 11},
      {6, 0, 7},
      {7, 0, 6},
      {2, 1, 9},
      {3, 1, 10},
      {6, 1, 8},
      {7, 1, 7},
      {2, 4, 1},
      {3, 4, 2},
      {6, 4, 10},
      {7, 4, 11},
      {2, 5, 0},
      {3, 5, 1},
      {6, 5, 9},
      {7, 5, 10}}},

    {4, /* (  0 ,  1 , -1 ) */
     {{2, 1, 11}, {6, 5, 11}, {2, 5, 2}, {6, 1, 6}}},

    {16, /* (  0 ,  0 ,  1 ) */
     {{1, 0, 12},
      {3, 0, 11},
      {5, 0, 5},
      {7, 0, 2},
      {1, 2, 9},
      {3, 2, 12},
      {5, 2, 8},
      {7, 2, 5},
      {1, 4, 3},
      {3, 4, 6},
      {5, 4, 12},
      {7, 4, 11},
      {1, 6, 0},
      {3, 6, 3},
      {5, 6, 9},
      {7, 6, 12}}}};

/**
 * @brief Get the size of the cell subtree.
 *
 * @param c The #cell.
 */
int cell_get_tree_size(struct cell *c) {
  /* Number of cells in this subtree. */
  int count = 1;

  /* Sum up the progeny if split. */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) count += cell_get_tree_size(c->progeny[k]);

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
  c->grav.parts_rebuild = gparts;

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
  c->stars.parts_rebuild = sparts;

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
 * @brief Link the cells recursively to the given #bpart array.
 *
 * @param c The #cell.
 * @param bparts The #bpart array.
 *
 * @return The number of particles linked.
 */
int cell_link_bparts(struct cell *c, struct bpart *bparts) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank)
    error("Linking foreign particles in a local cell!");

  if (c->black_holes.parts != NULL)
    error("Linking bparts into a cell that was already linked");
#endif

  c->black_holes.parts = bparts;

  /* Fill the progeny recursively, depth-first. */
  if (c->split) {
    int offset = 0;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        offset += cell_link_bparts(c->progeny[k], &bparts[offset]);
    }
  }

  /* Return the total number of linked particles. */
  return c->black_holes.count;
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
  if (cell_get_recv(c, task_subtype_xv) != NULL) {

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

  /* Do we have a gravity task at this level? */
  if (cell_get_recv(c, task_subtype_gpart) != NULL) {

    /* Recursively attach the gparts */
    const int counts = cell_link_gparts(c, gparts);
#ifdef SWIFT_DEBUG_CHECKS
    if (counts != c->grav.count)
      error("Something is wrong with the foreign counts");
#endif
    return counts;
  } else {
    c->grav.parts = gparts;
    c->grav.parts_rebuild = gparts;
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
 * @brief Recursively nullify all the particle pointers in a cell hierarchy.
 *
 * Should only be used on foreign cells!
 *
 * This will make any task or action running on these cells likely crash.
 * Recreating the foreign links will be necessary.
 *
 * @param c The #cell to act on.
 */
void cell_unlink_foreign_particles(struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank)
    error("Unlinking foreign particles in a local cell!");
#endif

  c->grav.parts = NULL;
  c->hydro.parts = NULL;
  c->stars.parts = NULL;
  c->black_holes.parts = NULL;

  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_unlink_foreign_particles(c->progeny[k]);
      }
    }
  }
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
  if (cell_get_recv(c, task_subtype_xv) != NULL) {
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

  /* Do we have a gravity task at this level? */
  if (cell_get_recv(c, task_subtype_gpart) != NULL) {
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
  float h_max_active = 0.f;
  float stars_h_max = 0.f;
  float stars_h_max_active = 0.f;

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
        h_max_active = max(h_max_active, c->progeny[k]->hydro.h_max_active);
        stars_h_max = max(stars_h_max, c->progeny[k]->stars.h_max);
        stars_h_max_active =
            max(stars_h_max_active, c->progeny[k]->stars.h_max_active);
      }
    }
  } else {
    /* Get the new value of h_max (note all particles are active) */
    for (int i = 0; i < count; ++i) h_max = max(h_max, parts[i].h);
    for (int i = 0; i < count; ++i)
      h_max_active = max(h_max_active, parts[i].h);
    for (int i = 0; i < scount; ++i)
      stars_h_max = max(stars_h_max, sparts[i].h);
    for (int i = 0; i < scount; ++i)
      stars_h_max_active = max(stars_h_max_active, sparts[i].h);
  }

  /* Record the change */
  c->hydro.h_max = h_max;
  c->hydro.h_max_active = h_max_active;
  c->stars.h_max = stars_h_max;
  c->stars.h_max_active = stars_h_max_active;
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
  c->rt.rt_gradient = NULL;
  c->rt.rt_transport = NULL;
  c->grav.grav = NULL;
  c->grav.mm = NULL;
  c->stars.density = NULL;
  c->stars.prepare1 = NULL;
  c->stars.prepare2 = NULL;
  c->stars.feedback = NULL;
  c->sinks.swallow = NULL;
  c->sinks.do_sink_swallow = NULL;
  c->sinks.do_gas_swallow = NULL;
  c->black_holes.density = NULL;
  c->black_holes.swallow = NULL;
  c->black_holes.do_gas_swallow = NULL;
  c->black_holes.do_bh_swallow = NULL;
  c->black_holes.feedback = NULL;
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
 * @brief Checks that the #sink in a cell are at the
 * current point in time
 *
 * Calls error() if the cell is not at the current time.
 *
 * @param c Cell to act upon
 * @param data The current time on the integer time-line
 */
void cell_check_sink_drift_point(struct cell *c, void *data) {
#ifdef SWIFT_DEBUG_CHECKS

  const integertime_t ti_drift = *(integertime_t *)data;

  /* Only check local cells */
  if (c->nodeID != engine_rank) return;

  /* Only check cells with content */
  if (c->sinks.count == 0) return;

  if (c->sinks.ti_old_part != ti_drift)
    error(
        "Cell in an incorrect time-zone! c->sinks.ti_old_part=%lld "
        "ti_drift=%lld",
        c->sinks.ti_old_part, ti_drift);

  for (int i = 0; i < c->sinks.count; ++i)
    if (c->sinks.parts[i].ti_drift != ti_drift &&
        c->sinks.parts[i].time_bin != time_bin_inhibited)
      error(
          "sink-part in an incorrect time-zone! sink->ti_drift=%lld "
          "ti_drift=%lld",
          c->sinks.parts[i].ti_drift, ti_drift);
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
 * @param grav_props The properties of the gravity scheme.
 */
void cell_make_multipoles(struct cell *c, integertime_t ti_current,
                          const struct gravity_props *const grav_props) {

  /* Reset everything */
  gravity_reset(c->grav.multipole);

  if (c->split) {

    /* Start by recursing */
    for (int k = 0; k < 8; ++k) {
      if (c->progeny[k] != NULL)
        cell_make_multipoles(c->progeny[k], ti_current, grav_props);
    }

    /* Compute CoM of all progenies */
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
    const double mass_inv = 1. / mass;
    c->grav.multipole->CoM[0] = CoM[0] * mass_inv;
    c->grav.multipole->CoM[1] = CoM[1] * mass_inv;
    c->grav.multipole->CoM[2] = CoM[2] * mass_inv;
    c->grav.multipole->m_pole.vel[0] = vel[0] * mass_inv;
    c->grav.multipole->m_pole.vel[1] = vel[1] * mass_inv;
    c->grav.multipole->m_pole.vel[2] = vel[2] * mass_inv;

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

    /* Compute the multipole power */
    gravity_multipole_compute_power(&c->grav.multipole->m_pole);

  } else {
    if (c->grav.count > 0) {

      gravity_P2M(c->grav.multipole, c->grav.parts, c->grav.count, grav_props);

      /* Compute the multipole power */
      gravity_multipole_compute_power(&c->grav.multipole->m_pole);

    } else {

      /* No gparts in that leaf cell */

      /* Set the values to something sensible */
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

    if (fabs(M_000 / c->grav.multipole->m_pole.M_000 - 1.) > 1e-2)
      error("Mass in progenies does not match!");
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
 * @param grav_props The properties of the gravity scheme.
 */
void cell_check_multipole(struct cell *c,
                          const struct gravity_props *const grav_props) {

#ifdef SWIFT_DEBUG_CHECKS
  struct gravity_tensors ma;
  const double tolerance = 1e-3; /* Relative */

  /* First recurse */
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        cell_check_multipole(c->progeny[k], grav_props);

  if (c->grav.count > 0) {
    /* Brute-force calculation */
    gravity_P2M(&ma, c->grav.parts, c->grav.count, grav_props);
    gravity_multipole_compute_power(&ma.m_pole);

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
  cell_clear_flag(c, cell_flag_do_hydro_drift | cell_flag_do_hydro_sub_drift |
                         cell_flag_do_grav_drift | cell_flag_do_grav_sub_drift |
                         cell_flag_do_bh_drift | cell_flag_do_bh_sub_drift |
                         cell_flag_do_stars_drift |
                         cell_flag_do_stars_sub_drift |
                         cell_flag_do_sink_drift | cell_flag_do_sink_sub_drift);
}

/**
 * @brief Clear the limiter flags on the given cell.
 */
void cell_clear_limiter_flags(struct cell *c, void *data) {
  cell_clear_flag(c,
                  cell_flag_do_hydro_limiter | cell_flag_do_hydro_sub_limiter);
}

/**
 * @brief Set the super-cell pointers for all cells in a hierarchy.
 *
 * @param c The top-level #cell to play with.
 * @param super Pointer to the deepest cell with tasks in this part of the
 * tree.
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
 * @param super_hydro Pointer to the deepest cell with tasks in this part of
 * the tree.
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
  return (c->timestep_collect != NULL || c->mpi.recv != NULL);
#else
  return (c->timestep_collect != NULL);
#endif
}

/**
 * @brief Resets all the sorting properties for the stars in a given cell
 * hierarchy.
 *
 * The clear_unused_flags argument can be used to additionally clean up all
 * the flags demanding a sort for the given cell. This should be used with
 * caution as it will prevent the sort tasks from doing anything on that cell
 * until these flags are reset.
 *
 * @param c The #cell to clean.
 * @param clear_unused_flags Do we also clean the flags demanding a sort?
 */
void cell_clear_stars_sort_flags(struct cell *c, const int clear_unused_flags) {

  /* Clear the flags that have not been reset by the sort task? */
  if (clear_unused_flags) {
    c->stars.requires_sorts = 0;
    c->stars.do_sort = 0;
    cell_clear_flag(c, cell_flag_do_stars_sub_sort);
  }

  /* Indicate that the cell is not sorted and cancel the pointer sorting
   * arrays.
   */
  c->stars.sorted = 0;
  cell_free_stars_sorts(c);

  /* Recurse if possible */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        cell_clear_stars_sort_flags(c->progeny[k], clear_unused_flags);
  }
}

/**
 * @brief Resets all the sorting properties for the hydro in a given cell
 * hierarchy.
 *
 * The clear_unused_flags argument can be used to additionally clean up all
 * the flags demanding a sort for the given cell. This should be used with
 * caution as it will prevent the sort tasks from doing anything on that cell
 * until these flags are reset.
 *
 * @param c The #cell to clean.
 * @param clear_unused_flags Do we also clean the flags demanding a sort?
 */
void cell_clear_hydro_sort_flags(struct cell *c, const int clear_unused_flags) {

  /* Clear the flags that have not been reset by the sort task? */
  if (clear_unused_flags) {
    c->hydro.do_sort = 0;
    c->hydro.requires_sorts = 0;
    cell_clear_flag(c, cell_flag_do_hydro_sub_sort);
  }

  /* Indicate that the cell is not sorted and cancel the pointer sorting
   * arrays.
   */
  c->hydro.sorted = 0;
  cell_free_hydro_sorts(c);

  /* Recurse if possible */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        cell_clear_hydro_sort_flags(c->progeny[k], clear_unused_flags);
  }
}

/**
 * @brief Recursively checks that all particles in a cell have a time-step
 */
void cell_check_timesteps(const struct cell *c, const integertime_t ti_current,
                          const timebin_t max_bin) {
#ifdef SWIFT_DEBUG_CHECKS

  if (c->hydro.ti_end_min == 0 && c->grav.ti_end_min == 0 &&
      c->stars.ti_end_min == 0 && c->black_holes.ti_end_min == 0 &&
      c->sinks.ti_end_min == 0 && c->rt.ti_rt_end_min == 0 && c->nr_tasks > 0)
    error("Cell without assigned time-step");

  if (c->split) {
    for (int k = 0; k < 8; ++k)
      if (c->progeny[k] != NULL)
        cell_check_timesteps(c->progeny[k], ti_current, max_bin);
  } else {
    if (c->nodeID == engine_rank)
      for (int i = 0; i < c->hydro.count; ++i)
        if (c->hydro.parts[i].time_bin == 0)
          error("Particle without assigned time-bin");
  }

  /* Other checks not relevent when starting-up */
  if (ti_current == 0) return;

  integertime_t ti_end_min = max_nr_timesteps;
  integertime_t ti_beg_max = 0;

  int count = 0;

  for (int i = 0; i < c->hydro.count; ++i) {

    const struct part *p = &c->hydro.parts[i];
    if (p->time_bin == time_bin_inhibited) continue;
    if (p->time_bin == time_bin_not_created) continue;

    ++count;

    integertime_t ti_end, ti_beg;

    if (p->time_bin <= max_bin) {
      integertime_t time_step = get_integer_timestep(p->time_bin);
      ti_end = get_integer_time_end(ti_current, p->time_bin) + time_step;
      ti_beg = get_integer_time_begin(ti_current + 1, p->time_bin);
    } else {
      ti_end = get_integer_time_end(ti_current, p->time_bin);
      ti_beg = get_integer_time_begin(ti_current + 1, p->time_bin);
    }

    ti_end_min = min(ti_end, ti_end_min);
    ti_beg_max = max(ti_beg, ti_beg_max);
  }

  /* Only check cells that have at least one non-inhibited particle */
  if (count > 0) {

    if (count != c->hydro.count) {

      /* Note that we use a < as the particle with the smallest time-bin
         might have been swallowed. This means we will run this cell with
         0 active particles but that's not wrong */
      if (ti_end_min < c->hydro.ti_end_min)
        error(
            "Non-matching ti_end_min. Cell=%lld true=%lld ti_current=%lld "
            "depth=%d",
            c->hydro.ti_end_min, ti_end_min, ti_current, c->depth);

    } else /* Normal case: nothing was swallowed/converted */ {
      if (ti_end_min != c->hydro.ti_end_min)
        error(
            "Non-matching ti_end_min. Cell=%lld true=%lld ti_current=%lld "
            "depth=%d",
            c->hydro.ti_end_min, ti_end_min, ti_current, c->depth);
    }

    if (ti_beg_max != c->hydro.ti_beg_max)
      error(
          "Non-matching ti_beg_max. Cell=%lld true=%lld ti_current=%lld "
          "depth=%d",
          c->hydro.ti_beg_max, ti_beg_max, ti_current, c->depth);
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
 * @brief Checks that a cell and all its progenies have cleared their sort
 * flags.
 *
 * Should only be used for debugging purposes.
 *
 * @param c The #cell to check.
 */
void cell_check_sort_flags(const struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  const int do_hydro_sub_sort = cell_get_flag(c, cell_flag_do_hydro_sub_sort);
  const int do_stars_sub_sort = cell_get_flag(c, cell_flag_do_stars_sub_sort);

  if (do_hydro_sub_sort)
    error(
        "cell %lld has a hydro sub_sort flag set. Node=%d depth=%d maxdepth=%d",
        c->cellID, c->nodeID, c->depth, c->maxdepth);

  if (do_stars_sub_sort)
    error(
        "cell %lld has a stars sub_sort flag set. Node=%d depth=%d maxdepth=%d",
        c->cellID, c->nodeID, c->depth, c->maxdepth);

  if (c->split) {
    for (int k = 0; k < 8; ++k) {
      if (c->progeny[k] != NULL) cell_check_sort_flags(c->progeny[k]);
    }
  }
#endif
}

/**
 * @brief Can we use the MM interactions fo a given pair of cells?
 *
 * The two cells have to be different!
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param e The #engine.
 * @param s The #space.
 * @param use_rebuild_data Are we considering the data at the last tree-build
 * (1) or current data (0)?
 * @param is_tree_walk Are we calling this in the tree walk (1) or for the
 * top-level task construction (0)?
 */
int cell_can_use_pair_mm(const struct cell *restrict ci,
                         const struct cell *restrict cj, const struct engine *e,
                         const struct space *s, const int use_rebuild_data,
                         const int is_tree_walk) {

  const struct gravity_props *props = e->gravity_properties;
  const int periodic = s->periodic;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};

  /* Check for trivial cases */
  if (is_tree_walk && ci->grav.count <= 1) return 0;
  if (is_tree_walk && cj->grav.count <= 1) return 0;

  /* Recover the multipole information */
  const struct gravity_tensors *restrict multi_i = ci->grav.multipole;
  const struct gravity_tensors *restrict multi_j = cj->grav.multipole;

  double dx, dy, dz;

  /* Get the distance between the CoMs */
  if (use_rebuild_data) {
    dx = multi_i->CoM_rebuild[0] - multi_j->CoM_rebuild[0];
    dy = multi_i->CoM_rebuild[1] - multi_j->CoM_rebuild[1];
    dz = multi_i->CoM_rebuild[2] - multi_j->CoM_rebuild[2];
  } else {
    dx = multi_i->CoM[0] - multi_j->CoM[0];
    dy = multi_i->CoM[1] - multi_j->CoM[1];
    dz = multi_i->CoM[2] - multi_j->CoM[2];
  }

  /* Apply BC */
  if (periodic) {
    dx = nearest(dx, dim[0]);
    dy = nearest(dy, dim[1]);
    dz = nearest(dz, dim[2]);
  }
  const double r2 = dx * dx + dy * dy + dz * dz;

  return gravity_M2L_accept_symmetric(props, multi_i, multi_j, r2,
                                      use_rebuild_data, periodic);
}
