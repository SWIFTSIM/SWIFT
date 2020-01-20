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
  c->hydro.rt_inject = NULL;
  c->grav.grav = NULL;
  c->grav.mm = NULL;
  c->stars.density = NULL;
  c->stars.feedback = NULL;
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
  if (c->timestep != NULL || c->mpi.recv != NULL) return 1;
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
  force = (force || cell_get_flag(c, cell_flag_do_hydro_drift));

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_part) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->hydro.count == 0) {
    /* Clear the drift flags. */
    cell_clear_flag(c, cell_flag_do_hydro_drift | cell_flag_do_hydro_sub_drift);

    /* Update the time of the last drift */
    c->hydro.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || cell_get_flag(c, cell_flag_do_hydro_sub_drift))) {

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

      /* Apply the effects of feedback on this particle
       * (Note: Only used in schemes that have a delayed feedback mechanism
       * otherwise just an empty function) */
      feedback_update_part(p, xp, e);

      /* Drift... */
      drift_part(p, xp, dt_drift, dt_kick_hydro, dt_kick_grav, dt_therm,
                 ti_old_part, ti_current, e->cosmology, e->hydro_properties,
                 e->entropy_floor);

      /* Update the tracers properties */
      tracers_after_drift(p, xp, e->internal_units, e->physical_constants,
                          with_cosmology, e->cosmology, e->hydro_properties,
                          e->cooling_func, e->time);

#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure the particle does not drift by more than a box length. */
      if (fabs(xp->v_full[0] * dt_drift) > e->s->dim[0] ||
          fabs(xp->v_full[1] * dt_drift) > e->s->dim[1] ||
          fabs(xp->v_full[2] * dt_drift) > e->s->dim[2]) {
        error(
            "Particle drifts by more than a box length! id %llu xp->v_full "
            "%.5e %.5e %.5e p->v %.5e %.5e %.5e",
            p->id, xp->v_full[0], xp->v_full[1], xp->v_full[2], p->v[0],
            p->v[1], p->v[2]);
      }
#endif

      /* In non-periodic BC runs, remove particles that crossed the border */
      if (!periodic) {

        /* Did the particle leave the box?  */
        if ((p->x[0] > dim[0]) || (p->x[0] < 0.) ||  // x
            (p->x[1] > dim[1]) || (p->x[1] < 0.) ||  // y
            (p->x[2] > dim[2]) || (p->x[2] < 0.)) {  // z

          lock_lock(&e->s->lock);

          /* Re-check that the particle has not been removed
           * by another thread before we do the deed. */
          if (!part_is_inhibited(p, e)) {

#ifdef WITH_LOGGER
            if (e->policy & engine_policy_logger) {
              /* Log the particle one last time. */
              logger_log_part(
                  e->logger, p, xp, e, /* log_all */ 1,
                  logger_pack_flags_and_data(logger_flag_delete, 0));
            }
#endif

            /* One last action before death? */
            hydro_remove_part(p, xp);

            /* Remove the particle entirely */
            cell_remove_part(e, c, p, xp);
          }

          if (lock_unlock(&e->s->lock) != 0)
            error("Failed to unlock the space!");

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

      /* Update the maximal smoothing length in the cell */
      cell_h_max = max(cell_h_max, p->h);

      /* Mark the particle has not being swallowed */
      black_holes_mark_part_as_not_swallowed(&p->black_holes_data);

      /* Get ready for a density calculation */
      if (part_is_active(p, e)) {
        hydro_init_part(p, &e->s->hs);
        black_holes_init_potential(&p->black_holes_data);
        chemistry_init_part(p, e->chemistry);
        pressure_floor_init_part(p, xp);
        star_formation_init_part(p, e->star_formation);
        tracers_after_init(p, xp, e->internal_units, e->physical_constants,
                           with_cosmology, e->cosmology, e->hydro_properties,
                           e->cooling_func, e->time);
        rt_init_part(p);

        /* Update the maximal active smoothing length in the cell */
        cell_h_max_active = max(cell_h_max_active, p->h);
      }
    }

    /* Now, get the maximal particle motion from its square */
    dx_max = sqrtf(dx2_max);
    dx_max_sort = sqrtf(dx2_max_sort);

    /* Store the values */
    c->hydro.h_max = cell_h_max;
    c->hydro.h_max_active = cell_h_max_active;
    c->hydro.dx_max_part = dx_max;
    c->hydro.dx_max_sort = dx_max_sort;

    /* Update the time of the last drift */
    c->hydro.ti_old_part = ti_current;
  }

  /* Clear the drift flags. */
  cell_clear_flag(c, cell_flag_do_hydro_drift | cell_flag_do_hydro_sub_drift);
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
  const struct gravity_props *grav_props = e->gravity_properties;

  /* Drift irrespective of cell flags? */
  force = (force || cell_get_flag(c, cell_flag_do_grav_drift));

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_gpart) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->grav.count == 0) {
    /* Clear the drift flags. */
    cell_clear_flag(c, cell_flag_do_grav_drift | cell_flag_do_grav_sub_drift);

    /* Update the time of the last drift */
    c->grav.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || cell_get_flag(c, cell_flag_do_grav_sub_drift))) {

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
      drift_gpart(gp, dt_drift, ti_old_gpart, ti_current, grav_props, e);

#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure the particle does not drift by more than a box length. */
      if (fabs(gp->v_full[0] * dt_drift) > e->s->dim[0] ||
          fabs(gp->v_full[1] * dt_drift) > e->s->dim[1] ||
          fabs(gp->v_full[2] * dt_drift) > e->s->dim[2]) {
        error(
            "Particle drifts by more than a box length! gp->v_full %.5e %.5e "
            "%.5e",
            gp->v_full[0], gp->v_full[1], gp->v_full[2]);
      }
#endif

      /* In non-periodic BC runs, remove particles that crossed the border */
      if (!periodic) {

        /* Did the particle leave the box?  */
        if ((gp->x[0] > dim[0]) || (gp->x[0] < 0.) ||  // x
            (gp->x[1] > dim[1]) || (gp->x[1] < 0.) ||  // y
            (gp->x[2] > dim[2]) || (gp->x[2] < 0.)) {  // z

          lock_lock(&e->s->lock);

          /* Re-check that the particle has not been removed
           * by another thread before we do the deed. */
          if (!gpart_is_inhibited(gp, e)) {

            /* Remove the particle entirely */
            if (gp->type == swift_type_dark_matter) {

#ifdef WITH_LOGGER
              if (e->policy & engine_policy_logger) {
                /* Log the particle one last time. */
                logger_log_gpart(
                    e->logger, gp, e, /* log_all */ 1,
                    logger_pack_flags_and_data(logger_flag_delete, 0));
              }
#endif

              /* Remove the particle */
              cell_remove_gpart(e, c, gp);
            }
          }

          if (lock_unlock(&e->s->lock) != 0)
            error("Failed to unlock the space!");

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
  cell_clear_flag(c, cell_flag_do_grav_drift | cell_flag_do_grav_sub_drift);
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
  force = (force || cell_get_flag(c, cell_flag_do_stars_drift));

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_spart) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->stars.count == 0) {
    /* Clear the drift flags. */
    cell_clear_flag(c, cell_flag_do_stars_drift | cell_flag_do_stars_sub_drift);

    /* Update the time of the last drift */
    c->stars.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || cell_get_flag(c, cell_flag_do_stars_sub_drift))) {

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

          lock_lock(&e->s->lock);

          /* Re-check that the particle has not been removed
           * by another thread before we do the deed. */
          if (!spart_is_inhibited(sp, e)) {

#ifdef WITH_LOGGER
            if (e->policy & engine_policy_logger) {
              /* Log the particle one last time. */
              logger_log_spart(
                  e->logger, sp, e, /* log_all */ 1,
                  logger_pack_flags_and_data(logger_flag_delete, 0));
            }
#endif

            /* Remove the particle entirely */
            cell_remove_spart(e, c, sp);
          }

          if (lock_unlock(&e->s->lock) != 0)
            error("Failed to unlock the space!");

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
        feedback_init_spart(sp);
        rt_init_spart(sp);
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
  cell_clear_flag(c, cell_flag_do_stars_drift | cell_flag_do_stars_sub_drift);
}

/**
 * @brief Recursively drifts the #bpart in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 * @param force Drift the particles irrespective of the #cell flags.
 */
void cell_drift_bpart(struct cell *c, const struct engine *e, int force) {

  const int periodic = e->s->periodic;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const float black_holes_h_max = e->hydro_properties->h_max;
  const float black_holes_h_min = e->hydro_properties->h_min;
  const integertime_t ti_old_bpart = c->black_holes.ti_old_part;
  const integertime_t ti_current = e->ti_current;
  struct bpart *const bparts = c->black_holes.parts;

  float dx_max = 0.f, dx2_max = 0.f;
  float cell_h_max = 0.f;

  /* Drift irrespective of cell flags? */
  force = (force || cell_get_flag(c, cell_flag_do_bh_drift));

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_bpart) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->black_holes.count == 0) {

    /* Clear the drift flags. */
    cell_clear_flag(c, cell_flag_do_bh_drift | cell_flag_do_bh_sub_drift);

    /* Update the time of the last drift */
    c->black_holes.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || cell_get_flag(c, cell_flag_do_bh_sub_drift))) {

    /* Loop over the progeny and collect their data. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];

        /* Recurse */
        cell_drift_bpart(cp, e, force);

        /* Update */
        dx_max = max(dx_max, cp->black_holes.dx_max_part);
        cell_h_max = max(cell_h_max, cp->black_holes.h_max);
      }
    }

    /* Store the values */
    c->black_holes.h_max = cell_h_max;
    c->black_holes.dx_max_part = dx_max;

    /* Update the time of the last drift */
    c->black_holes.ti_old_part = ti_current;

  } else if (!c->split && force && ti_current > ti_old_bpart) {

    /* Drift from the last time the cell was drifted to the current time */
    double dt_drift;
    if (with_cosmology) {
      dt_drift =
          cosmology_get_drift_factor(e->cosmology, ti_old_bpart, ti_current);
    } else {
      dt_drift = (ti_current - ti_old_bpart) * e->time_base;
    }

    /* Loop over all the star particles in the cell */
    const size_t nr_bparts = c->black_holes.count;
    for (size_t k = 0; k < nr_bparts; k++) {

      /* Get a handle on the bpart. */
      struct bpart *const bp = &bparts[k];

      /* Ignore inhibited particles */
      if (bpart_is_inhibited(bp, e)) continue;

      /* Drift... */
      drift_bpart(bp, dt_drift, ti_old_bpart, ti_current);

#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure the particle does not drift by more than a box length. */
      if (fabs(bp->v[0] * dt_drift) > e->s->dim[0] ||
          fabs(bp->v[1] * dt_drift) > e->s->dim[1] ||
          fabs(bp->v[2] * dt_drift) > e->s->dim[2]) {
        error("Particle drifts by more than a box length!");
      }
#endif

      /* In non-periodic BC runs, remove particles that crossed the border */
      if (!periodic) {

        /* Did the particle leave the box?  */
        if ((bp->x[0] > dim[0]) || (bp->x[0] < 0.) ||  // x
            (bp->x[1] > dim[1]) || (bp->x[1] < 0.) ||  // y
            (bp->x[2] > dim[2]) || (bp->x[2] < 0.)) {  // z

          lock_lock(&e->s->lock);

          /* Re-check that the particle has not been removed
           * by another thread before we do the deed. */
          if (!bpart_is_inhibited(bp, e)) {

#ifdef WITH_LOGGER
            if (e->policy & engine_policy_logger) {
              error("Logging of black hole particles is not yet implemented.");
            }
#endif

            /* Remove the particle entirely */
            cell_remove_bpart(e, c, bp);
          }

          if (lock_unlock(&e->s->lock) != 0)
            error("Failed to unlock the space!");

          continue;
        }
      }

      /* Limit h to within the allowed range */
      bp->h = min(bp->h, black_holes_h_max);
      bp->h = max(bp->h, black_holes_h_min);

      /* Compute (square of) motion since last cell construction */
      const float dx2 = bp->x_diff[0] * bp->x_diff[0] +
                        bp->x_diff[1] * bp->x_diff[1] +
                        bp->x_diff[2] * bp->x_diff[2];
      dx2_max = max(dx2_max, dx2);

      /* Maximal smoothing length */
      cell_h_max = max(cell_h_max, bp->h);

      /* Mark the particle has not being swallowed */
      black_holes_mark_bpart_as_not_swallowed(&bp->merger_data);

      /* Get ready for a density calculation */
      if (bpart_is_active(bp, e)) {
        black_holes_init_bpart(bp);
      }
    }

    /* Now, get the maximal particle motion from its square */
    dx_max = sqrtf(dx2_max);

    /* Store the values */
    c->black_holes.h_max = cell_h_max;
    c->black_holes.dx_max_part = dx_max;

    /* Update the time of the last drift */
    c->black_holes.ti_old_part = ti_current;
  }

  /* Clear the drift flags. */
  cell_clear_flag(c, cell_flag_do_bh_drift | cell_flag_do_bh_sub_drift);
}

/**
 * @brief Recursively drifts the #sink's in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 * @param force Drift the particles irrespective of the #cell flags.
 */
void cell_drift_sink(struct cell *c, const struct engine *e, int force) {

  const int periodic = e->s->periodic;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const integertime_t ti_old_sink = c->sinks.ti_old_part;
  const integertime_t ti_current = e->ti_current;
  struct sink *const sinks = c->sinks.parts;

  float dx_max = 0.f, dx2_max = 0.f;
  float cell_r_max = 0.f;

  /* Drift irrespective of cell flags? */
  force = (force || cell_get_flag(c, cell_flag_do_sink_drift));

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_sink) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->sinks.count == 0) {

    /* Clear the drift flags. */
    cell_clear_flag(c, cell_flag_do_sink_drift | cell_flag_do_sink_sub_drift);

    /* Update the time of the last drift */
    c->sinks.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || cell_get_flag(c, cell_flag_do_sink_sub_drift))) {

    /* Loop over the progeny and collect their data. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];

        /* Recurse */
        cell_drift_sink(cp, e, force);

        /* Update */
        dx_max = max(dx_max, cp->sinks.dx_max_part);
        cell_r_max = max(cell_r_max, cp->sinks.r_cut_max);
      }
    }

    /* Store the values */
    c->sinks.r_cut_max = cell_r_max;
    c->sinks.dx_max_part = dx_max;

    /* Update the time of the last drift */
    c->sinks.ti_old_part = ti_current;

  } else if (!c->split && force && ti_current > ti_old_sink) {

    /* Drift from the last time the cell was drifted to the current time */
    double dt_drift;
    if (with_cosmology) {
      dt_drift =
          cosmology_get_drift_factor(e->cosmology, ti_old_sink, ti_current);
    } else {
      dt_drift = (ti_current - ti_old_sink) * e->time_base;
    }

    /* Loop over all the star particles in the cell */
    const size_t nr_sinks = c->sinks.count;
    for (size_t k = 0; k < nr_sinks; k++) {

      /* Get a handle on the sink. */
      struct sink *const sink = &sinks[k];

      /* Ignore inhibited particles */
      if (sink_is_inhibited(sink, e)) continue;

      /* Drift... */
      drift_sink(sink, dt_drift, ti_old_sink, ti_current);

#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure the particle does not drift by more than a box length. */
      if (fabs(sink->v[0] * dt_drift) > e->s->dim[0] ||
          fabs(sink->v[1] * dt_drift) > e->s->dim[1] ||
          fabs(sink->v[2] * dt_drift) > e->s->dim[2]) {
        error("Particle drifts by more than a box length!");
      }
#endif

      /* In non-periodic BC runs, remove particles that crossed the border */
      if (!periodic) {

        /* Did the particle leave the box?  */
        if ((sink->x[0] > dim[0]) || (sink->x[0] < 0.) ||  // x
            (sink->x[1] > dim[1]) || (sink->x[1] < 0.) ||  // y
            (sink->x[2] > dim[2]) || (sink->x[2] < 0.)) {  // z

          lock_lock(&e->s->lock);

          /* Re-check that the particle has not been removed
           * by another thread before we do the deed. */
          if (!sink_is_inhibited(sink, e)) {

#ifdef WITH_LOGGER
            if (e->policy & engine_policy_logger) {
              error("Logging of sink particles is not yet implemented.");
            }
#endif

            /* Remove the particle entirely */
            // cell_remove_sink(e, c, bp);
            error("TODO: loic implement cell_remove_sink");
          }

          if (lock_unlock(&e->s->lock) != 0)
            error("Failed to unlock the space!");

          continue;
        }
      }

      /* sp->h does not need to be limited. */

      /* Compute (square of) motion since last cell construction */
      const float dx2 = sink->x_diff[0] * sink->x_diff[0] +
                        sink->x_diff[1] * sink->x_diff[1] +
                        sink->x_diff[2] * sink->x_diff[2];
      dx2_max = max(dx2_max, dx2);

      /* Maximal smoothing length */
      cell_r_max = max(cell_r_max, sink->r_cut);

      /* Get ready for a density calculation */
      if (sink_is_active(sink, e)) {
        sink_init_sink(sink);
      }
    }

    /* Now, get the maximal particle motion from its square */
    dx_max = sqrtf(dx2_max);

    /* Store the values */
    c->sinks.r_cut_max = cell_r_max;
    c->sinks.dx_max_part = dx_max;

    /* Update the time of the last drift */
    c->sinks.ti_old_part = ti_current;
  }

  /* Clear the drift flags. */
  cell_clear_flag(c, cell_flag_do_sink_drift | cell_flag_do_sink_sub_drift);
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
      c->sinks.ti_end_min == 0 && c->nr_tasks > 0)
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
  integertime_t ti_end_max = 0;
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
    ti_end_max = max(ti_end, ti_end_max);
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

    if (ti_end_max > c->hydro.ti_end_max)
      error(
          "Non-matching ti_end_max. Cell=%lld true=%lld ti_current=%lld "
          "depth=%d",
          c->hydro.ti_end_max, ti_end_max, ti_current, c->depth);

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
