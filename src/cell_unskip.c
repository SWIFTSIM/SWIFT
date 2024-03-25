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
#include "active.h"
#include "engine.h"
#include "feedback.h"
#include "space_getsid.h"

extern int engine_star_resort_task_depth;

/**
 * @brief Recursively clear the stars_resort flag in a cell hierarchy.
 *
 * @param c The #cell to act on.
 */
void cell_set_star_resort_flag(struct cell *c) {

  cell_set_flag(c, cell_flag_do_stars_resort);

  /* Abort if we reched the level where the resorting task lives */
  if (c->depth == engine_star_resort_task_depth || c->hydro.super == c) return;

  if (c->split) {
    for (int k = 0; k < 8; ++k)
      if (c->progeny[k] != NULL) cell_set_star_resort_flag(c->progeny[k]);
  }
}

/**
 * @brief Recurses in a cell hierarchy down to the level where the
 * star resort tasks are and activates them.
 *
 * The function will fail if called *below* the super-level
 *
 * @param c The #cell to recurse into.
 * @param s The #scheduler.
 */
void cell_activate_star_resort_tasks(struct cell *c, struct scheduler *s) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->hydro.super != NULL && c->hydro.super != c)
    error("Function called below the super level!");
#endif

  /* The resort tasks are at either the chosen depth or the super level,
   * whichever comes first. */
  if ((c->depth == engine_star_resort_task_depth || c->hydro.super == c) &&
      c->hydro.count > 0) {
    scheduler_activate(s, c->hydro.stars_resort);
  } else {
    for (int k = 0; k < 8; ++k) {
      if (c->progeny[k] != NULL) {
        cell_activate_star_resort_tasks(c->progeny[k], s);
      }
    }
  }
}

/**
 * @brief Activate the star formation task as well as the resorting of stars
 *
 * Must be called at the top-level in the tree (where the SF task is...)
 *
 * @param c The (top-level) #cell.
 * @param s The #scheduler.
 * @param with_feedback Are we running with feedback?
 */
void cell_activate_star_formation_tasks(struct cell *c, struct scheduler *s,
                                        const int with_feedback) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->depth != 0) error("Function should be called at the top-level only");
#endif

  /* Have we already unskipped that task? */
  if (c->hydro.star_formation->skip == 0) return;

  /* Activate the star formation task */
  scheduler_activate(s, c->hydro.star_formation);

  /* Activate the star resort tasks at whatever level they are */
  if (with_feedback) {
    cell_activate_star_resort_tasks(c, s);
  }
}

/**
 * @brief Activate the star formation task from the sink as well as the
 * resorting of stars
 *
 * Must be called at the top-level in the tree (where the SF task is...)
 *
 * @param c The (top-level) #cell.
 * @param s The #scheduler.
 * @param with_feedback Are we running with feedback?
 */
void cell_activate_star_formation_sink_tasks(struct cell *c,
                                             struct scheduler *s,
                                             const int with_feedback) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->depth != 0) error("Function should be called at the top-level only");
#endif

  /* Have we already unskipped that task? */
  if (c->sinks.star_formation_sink->skip == 0) return;

  /* Activate the star formation task */
  scheduler_activate(s, c->sinks.star_formation_sink);

  /* Activate the star resort tasks at whatever level they are */
  if (with_feedback) {
    cell_activate_star_resort_tasks(c, s);
  }
}

/**
 * @brief Activate the sink formation task.
 *
 * Must be called at the top-level in the tree (where the SF task is...)
 *
 * @param c The (top-level) #cell.
 * @param s The #scheduler.
 */
void cell_activate_sink_formation_tasks(struct cell *c, struct scheduler *s) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->depth != 0) error("Function should be called at the top-level only");
#endif

  /* Have we already unskipped that task? */
  if (c->sinks.sink_formation->skip == 0) return;

  /* Activate the star formation task */
  scheduler_activate(s, c->sinks.sink_formation);
}

/**
 * @brief Recursively activate the hydro ghosts (and implicit links) in a cell
 * hierarchy.
 *
 * @param c The #cell.
 * @param s The #scheduler.
 * @param e The #engine.
 */
void cell_recursively_activate_hydro_ghosts(struct cell *c, struct scheduler *s,
                                            const struct engine *e) {
  /* Early abort? */
  if ((c->hydro.count == 0) || !cell_is_active_hydro(c, e)) return;

  /* Is the ghost at this level? */
  if (c->hydro.ghost != NULL) {
    scheduler_activate(s, c->hydro.ghost);
  } else {

#ifdef SWIFT_DEBUG_CHECKS
    if (!c->split)
      error("Reached the leaf level without finding a hydro ghost!");
#endif

    /* Keep recursing */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        cell_recursively_activate_hydro_ghosts(c->progeny[k], s, e);
  }
}

/**
 * @brief Activate the hydro ghosts (and implicit links) in a cell hierarchy.
 *
 * @param c The #cell.
 * @param s The #scheduler.
 * @param e The #engine.
 */
void cell_activate_hydro_ghosts(struct cell *c, struct scheduler *s,
                                const struct engine *e) {
  scheduler_activate(s, c->hydro.ghost_in);
  scheduler_activate(s, c->hydro.ghost_out);
  cell_recursively_activate_hydro_ghosts(c, s, e);
}

/**
 * @brief Recursively activate the cooling (and implicit links) in a cell
 * hierarchy.
 *
 * @param c The #cell.
 * @param s The #scheduler.
 * @param e The #engine.
 */
void cell_recursively_activate_cooling(struct cell *c, struct scheduler *s,
                                       const struct engine *e) {
  /* Early abort? */
  if ((c->hydro.count == 0) || !cell_is_active_hydro(c, e)) return;

  /* Is the ghost at this level? */
  if (c->hydro.cooling != NULL) {
    scheduler_activate(s, c->hydro.cooling);
  } else {

#ifdef SWIFT_DEBUG_CHECKS
    if (!c->split)
      error("Reached the leaf level without finding a cooling task!");
#endif

    /* Keep recursing */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        cell_recursively_activate_cooling(c->progeny[k], s, e);
  }
}

/**
 * @brief Activate the cooling tasks (and implicit links) in a cell hierarchy.
 *
 * @param c The #cell.
 * @param s The #scheduler.
 * @param e The #engine.
 */
void cell_activate_cooling(struct cell *c, struct scheduler *s,
                           const struct engine *e) {
  scheduler_activate(s, c->hydro.cooling_in);
  scheduler_activate(s, c->hydro.cooling_out);
  cell_recursively_activate_cooling(c, s, e);
}

/**
 * @brief Recurse down in a cell hierarchy until the hydro.super level is
 * reached and activate the spart drift at that level.
 *
 * @param c The #cell to recurse into.
 * @param s The #scheduler.
 */
void cell_activate_super_spart_drifts(struct cell *c, struct scheduler *s) {

  /* Early abort?
   * We can stop if there is no gas as none of it will turn into stars */
  if (c->hydro.count == 0) return;

  if (c == c->hydro.super) {
    cell_activate_drift_spart(c, s);
    cell_set_flag(c, cell_flag_do_stars_drift);
  } else {
    if (c->split) {
      for (int k = 0; k < 8; ++k) {
        if (c->progeny[k] != NULL) {
          cell_activate_super_spart_drifts(c->progeny[k], s);
        }
      }
    } else {
#ifdef SWIFT_DEBUG_CHECKS
      error("Reached a leaf cell without finding a hydro.super!!");
#endif
    }
  }
}

/**
 * @brief Recurse down in a cell hierarchy until the hydro.super level is
 * reached and activate the sink drift at that level.
 *
 * @param c The #cell to recurse into.
 * @param s The #scheduler.
 */
void cell_activate_super_sink_drifts(struct cell *c, struct scheduler *s) {

  /* Early abort? */
  if (c->hydro.count == 0) return;

  if (c == c->hydro.super) {
    cell_activate_drift_sink(c, s);
  } else {
    if (c->split) {
      for (int k = 0; k < 8; ++k) {
        if (c->progeny[k] != NULL) {
          cell_activate_super_sink_drifts(c->progeny[k], s);
        }
      }
    } else {
#ifdef SWIFT_DEBUG_CHECKS
      error("Reached a leaf cell without finding a hydro.super!!");
#endif
    }
  }
}

/**
 * @brief Activate the #part drifts on the given cell.
 */
void cell_activate_drift_part(struct cell *c, struct scheduler *s) {
  /* If this cell is already marked for drift, quit early. */
  if (cell_get_flag(c, cell_flag_do_hydro_drift)) return;

  /* Mark this cell for drifting. */
  cell_set_flag(c, cell_flag_do_hydro_drift);

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
         parent != NULL && !cell_get_flag(parent, cell_flag_do_hydro_sub_drift);
         parent = parent->parent) {
      /* Mark this cell for drifting */
      cell_set_flag(parent, cell_flag_do_hydro_sub_drift);

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
 * @brief Activate the #part sync tasks on the given cell.
 */
void cell_activate_sync_part(struct cell *c, struct scheduler *s) {

  /* If this cell is already marked for sync, quit early. */
  if (cell_get_flag(c, cell_flag_do_hydro_sync)) return;

  /* Mark this cell for synchronization. */
  cell_set_flag(c, cell_flag_do_hydro_sync);

  /* Set the do_sub_sync all the way up and activate the super sync
     if this has not yet been done. */
  if (c == c->super) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->timestep_sync == NULL)
      error("Trying to activate un-existing c->timestep_sync");
#endif
    scheduler_activate(s, c->timestep_sync);
    scheduler_activate(s, c->top->timestep_collect);
    scheduler_activate(s, c->kick1);
  } else {
    for (struct cell *parent = c->parent;
         parent != NULL && !cell_get_flag(parent, cell_flag_do_hydro_sub_sync);
         parent = parent->parent) {
      /* Mark this cell for drifting */
      cell_set_flag(parent, cell_flag_do_hydro_sub_sync);

      if (parent == c->super) {
#ifdef SWIFT_DEBUG_CHECKS
        if (parent->timestep_sync == NULL)
          error("Trying to activate un-existing parent->timestep_sync");
#endif
        scheduler_activate(s, parent->timestep_sync);
        scheduler_activate(s, parent->top->timestep_collect);
        scheduler_activate(s, parent->kick1);
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
  if (cell_get_flag(c, cell_flag_do_grav_drift)) return;

  /* Mark this cell for drifting. */
  cell_set_flag(c, cell_flag_do_grav_drift);

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
         parent != NULL && !cell_get_flag(parent, cell_flag_do_grav_sub_drift);
         parent = parent->parent) {
      cell_set_flag(parent, cell_flag_do_grav_sub_drift);

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
  if (cell_get_flag(c, cell_flag_do_stars_drift)) return;

  /* Mark this cell for drifting. */
  cell_set_flag(c, cell_flag_do_stars_drift);

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
         parent != NULL && !cell_get_flag(parent, cell_flag_do_stars_sub_drift);
         parent = parent->parent) {
      /* Mark this cell for drifting */
      cell_set_flag(parent, cell_flag_do_stars_sub_drift);

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
 * @brief Activate the #bpart drifts on the given cell.
 */
void cell_activate_drift_bpart(struct cell *c, struct scheduler *s) {

  /* If this cell is already marked for drift, quit early. */
  if (cell_get_flag(c, cell_flag_do_bh_drift)) return;

  /* Mark this cell for drifting. */
  cell_set_flag(c, cell_flag_do_bh_drift);

  /* Set the do_black_holes_sub_drifts all the way up and activate the super
     drift if this has not yet been done. */
  if (c == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->black_holes.drift == NULL)
      error("Trying to activate un-existing c->black_holes.drift");
#endif
    scheduler_activate(s, c->black_holes.drift);
  } else {
    for (struct cell *parent = c->parent;
         parent != NULL && !cell_get_flag(parent, cell_flag_do_bh_sub_drift);
         parent = parent->parent) {
      /* Mark this cell for drifting */
      cell_set_flag(parent, cell_flag_do_bh_sub_drift);

      if (parent == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
        if (parent->black_holes.drift == NULL)
          error("Trying to activate un-existing parent->black_holes.drift");
#endif
        scheduler_activate(s, parent->black_holes.drift);
        break;
      }
    }
  }
}

/**
 * @brief Activate the #sink drifts on the given cell.
 */
void cell_activate_drift_sink(struct cell *c, struct scheduler *s) {

  /* If this cell is already marked for drift, quit early. */
  if (cell_get_flag(c, cell_flag_do_sink_drift)) return;

  /* Mark this cell for drifting. */
  cell_set_flag(c, cell_flag_do_sink_drift);

  /* Set the do_sink_sub_drifts all the way up and activate the super
     drift if this has not yet been done. */
  if (c == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->sinks.drift == NULL)
      error("Trying to activate un-existing c->sinks.drift");
#endif
    scheduler_activate(s, c->sinks.drift);
  } else {
    for (struct cell *parent = c->parent;
         parent != NULL && !cell_get_flag(parent, cell_flag_do_sink_sub_drift);
         parent = parent->parent) {
      /* Mark this cell for drifting */
      cell_set_flag(parent, cell_flag_do_sink_sub_drift);

      if (parent == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
        if (parent->sinks.drift == NULL)
          error("Trying to activate un-existing parent->sinks.drift");
#endif
        scheduler_activate(s, parent->sinks.drift);
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
  if (cell_get_flag(c, cell_flag_do_hydro_limiter)) return;

  /* Mark this cell for limiting. */
  cell_set_flag(c, cell_flag_do_hydro_limiter);

  /* Set the do_sub_limiter all the way up and activate the super limiter
     if this has not yet been done. */
  if (c == c->super) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->timestep_limiter == NULL)
      error("Trying to activate un-existing c->timestep_limiter");
#endif
    scheduler_activate(s, c->timestep_limiter);
    scheduler_activate(s, c->top->timestep_collect);
    scheduler_activate(s, c->kick1);
  } else {
    for (struct cell *parent = c->parent;
         parent != NULL &&
         !cell_get_flag(parent, cell_flag_do_hydro_sub_limiter);
         parent = parent->parent) {
      /* Mark this cell for limiting */
      cell_set_flag(parent, cell_flag_do_hydro_sub_limiter);

      if (parent == c->super) {
#ifdef SWIFT_DEBUG_CHECKS
        if (parent->timestep_limiter == NULL)
          error("Trying to activate un-existing parent->timestep_limiter");
#endif
        scheduler_activate(s, parent->timestep_limiter);
        scheduler_activate(s, parent->top->timestep_collect);
        scheduler_activate(s, parent->kick1);
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
    cell_set_flag(c, cell_flag_skip_rt_sort);
    if (c->nodeID == engine_rank) cell_activate_drift_part(c, s);
  } else {
    for (struct cell *parent = c->parent;
         parent != NULL && !cell_get_flag(parent, cell_flag_do_hydro_sub_sort);
         parent = parent->parent) {
      cell_set_flag(parent, cell_flag_do_hydro_sub_sort);
      if (parent == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
        if (parent->hydro.sorts == NULL)
          error("Trying to activate un-existing parents->hydro.sorts");
#endif
        scheduler_activate(s, parent->hydro.sorts);
        cell_set_flag(parent, cell_flag_skip_rt_sort);
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
 * @brief Activate the sorts up a cell hierarchy. Activate drifts
 * and hydro sorts on local cells, and rt_sorts on foreign cells.
 */
void cell_activate_rt_sorts_up(struct cell *c, struct scheduler *s) {

  cell_set_flag(c, cell_flag_rt_requests_sort);

  if (c == c->hydro.super) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->nodeID == engine_rank && c->hydro.sorts == NULL)
      error("Trying to activate non-existing c->hydro.sorts");
    if (c->nodeID != engine_rank && c->rt.rt_sorts == NULL)
      error("Trying to activate non-existing c->rt.rt_sorts");
#endif
    if (c->nodeID == engine_rank) {
      cell_set_flag(c, cell_flag_skip_rt_sort);
      scheduler_activate(s, c->hydro.sorts);
    } else {
      scheduler_activate(s, c->rt.rt_sorts);
    }
  } else {

    for (struct cell *parent = c->parent;
         parent != NULL && !cell_get_flag(parent, cell_flag_do_rt_sub_sort);
         parent = parent->parent) {

      /* Need a separate flag for RT sub sorts because RT sorts don't
       * necessarily activate the hydro sorts tasks, yet the do_hydro_sub_sort
       * flag is used as an early exit while climbing up the tree. */
      cell_set_flag(parent, cell_flag_do_rt_sub_sort);
      cell_set_flag(parent, cell_flag_rt_requests_sort);

      if (parent == c->hydro.super) {

#ifdef SWIFT_DEBUG_CHECKS
        if (parent->nodeID == engine_rank && parent->hydro.sorts == NULL)
          error("Trying to activate non-existing parents->hydro.sorts");
        if (parent->nodeID != engine_rank && parent->rt.rt_sorts == NULL)
          error("Trying to activate non-existing parents->rt.rt_sorts");
#endif

        if (parent->nodeID == engine_rank) {
          /* Mark the progeny to skip the RT sort as well */
          cell_set_flag(c, cell_flag_skip_rt_sort);
          scheduler_activate(s, parent->hydro.sorts);
        } else {
          scheduler_activate(s, parent->rt.rt_sorts);
        }
        break;
      }
    }
  }
}

/**
 * @brief Activate the sorts on a given cell, if needed. Activate
 * hydro sorts on local cells, and rt_sorts on foreign cells.
 */
void cell_activate_rt_sorts(struct cell *c, int sid, struct scheduler *s) {
  /* Do we need to re-sort? */
  if (c->hydro.dx_max_sort > space_maxreldx * c->dmin) {
    /* Climb up the tree to active the sorts in that direction */
    for (struct cell *finger = c; finger != NULL; finger = finger->parent) {
      if (finger->hydro.requires_sorts) {
        atomic_or(&finger->hydro.do_sort, finger->hydro.requires_sorts);
        cell_activate_rt_sorts_up(finger, s);
      }
      finger->hydro.sorted = 0;
    }
  }

  /* Has this cell been sorted at all for the given sid? */
  if (!(c->hydro.sorted & (1 << sid)) || c->nodeID != engine_rank) {
    atomic_or(&c->hydro.do_sort, (1 << sid));
    cell_activate_rt_sorts_up(c, s);
  }
}

/**
 * @brief Mark cells up a hierarchy to not run RT sorts.
 * */
void cell_set_skip_rt_sort_flag_up(struct cell *c) {

  for (struct cell *finger = c; finger != NULL; finger = finger->parent) {
    cell_set_flag(finger, cell_flag_skip_rt_sort);
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
      cell_activate_drift_spart(c, s);
    }
  } else {

    /* Climb up the tree and set the flags */
    for (struct cell *parent = c->parent;
         parent != NULL && !cell_get_flag(parent, cell_flag_do_stars_sub_sort);
         parent = parent->parent) {

      cell_set_flag(parent, cell_flag_do_stars_sub_sort);

      /* Reached the super-level? Activate the task and abort */
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
 * @param with_timestep_limiter Are we running with time-step limiting on?
 */
void cell_activate_subcell_hydro_tasks(struct cell *ci, struct cell *cj,
                                       struct scheduler *s,
                                       const int with_timestep_limiter) {
  const struct engine *e = s->space->e;

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
          cell_activate_subcell_hydro_tasks(ci->progeny[j], NULL, s,
                                            with_timestep_limiter);
          for (int k = j + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              cell_activate_subcell_hydro_tasks(ci->progeny[j], ci->progeny[k],
                                                s, with_timestep_limiter);
        }
      }
    } else {
      /* We have reached the bottom of the tree: activate drift */
      cell_activate_drift_part(ci, s);
      if (with_timestep_limiter) cell_activate_limiter(ci, s);
    }
  }

  /* Otherwise, pair interation */
  else {
    /* Should we even bother? */
    if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;
    if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

    /* Get the orientation of the pair. */
    double shift[3];
    const int sid = space_getsid(s->space, &ci, &cj, shift);

    /* recurse? */
    if (cell_can_recurse_in_pair_hydro_task(ci) &&
        cell_can_recurse_in_pair_hydro_task(cj)) {
      const struct cell_split_pair *csp = &cell_split_pairs[sid];
      for (int k = 0; k < csp->count; k++) {
        const int pid = csp->pairs[k].pid;
        const int pjd = csp->pairs[k].pjd;
        if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
          cell_activate_subcell_hydro_tasks(ci->progeny[pid], cj->progeny[pjd],
                                            s, with_timestep_limiter);
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
      if (ci->nodeID == engine_rank && with_timestep_limiter)
        cell_activate_limiter(ci, s);
      if (cj->nodeID == engine_rank && with_timestep_limiter)
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
 * @param with_star_formation Are we running with star formation switched on?
 * @param with_star_formation Are we running with star formation from the sinks
 * switched on?
 * @param with_timestep_sync Are we running with time-step synchronization on?
 */
void cell_activate_subcell_stars_tasks(struct cell *ci, struct cell *cj,
                                       struct scheduler *s,
                                       const int with_star_formation,
                                       const int with_star_formation_sink,
                                       const int with_timestep_sync) {
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

    const int ci_active = cell_need_activating_stars(ci, e, with_star_formation,
                                                     with_star_formation_sink);

    /* Do anything? */
    if (!ci_active || ci->hydro.count == 0 ||
        (!with_star_formation && !with_star_formation_sink &&
         ci->stars.count == 0))
      return;

    /* Recurse? */
    if (cell_can_recurse_in_self_stars_task(ci)) {
      /* Loop over all progenies and pairs of progenies */
      for (int j = 0; j < 8; j++) {
        if (ci->progeny[j] != NULL) {
          cell_activate_subcell_stars_tasks(
              ci->progeny[j], NULL, s, with_star_formation,
              with_star_formation_sink, with_timestep_sync);
          for (int k = j + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              cell_activate_subcell_stars_tasks(
                  ci->progeny[j], ci->progeny[k], s, with_star_formation,
                  with_star_formation_sink, with_timestep_sync);
        }
      }
    } else {
      /* We have reached the bottom of the tree: activate drift */
      cell_activate_drift_spart(ci, s);
      cell_activate_drift_part(ci, s);
      if (with_timestep_sync) cell_activate_sync_part(ci, s);
    }
  }

  /* Otherwise, pair interaction */
  else {

    /* Get the orientation of the pair. */
    double shift[3];
    const int sid = space_getsid(s->space, &ci, &cj, shift);

    const int ci_active = cell_need_activating_stars(ci, e, with_star_formation,
                                                     with_star_formation_sink);
    const int cj_active = cell_need_activating_stars(cj, e, with_star_formation,
                                                     with_star_formation_sink);

    /* Should we even bother? */
    if (!ci_active && !cj_active) return;

    /* recurse? */
    if (cell_can_recurse_in_pair_stars_task(ci, cj) &&
        cell_can_recurse_in_pair_stars_task(cj, ci)) {

      const struct cell_split_pair *csp = &cell_split_pairs[sid];
      for (int k = 0; k < csp->count; k++) {
        const int pid = csp->pairs[k].pid;
        const int pjd = csp->pairs[k].pjd;
        if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
          cell_activate_subcell_stars_tasks(
              ci->progeny[pid], cj->progeny[pjd], s, with_star_formation,
              with_star_formation_sink, with_timestep_sync);
      }
    }

    /* Otherwise, activate the sorts and drifts. */
    else {

      if (ci_active) {

        /* We are going to interact this pair, so store some values. */
        atomic_or(&cj->hydro.requires_sorts, 1 << sid);
        atomic_or(&ci->stars.requires_sorts, 1 << sid);

        cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;
        ci->stars.dx_max_sort_old = ci->stars.dx_max_sort;

        /* Activate the drifts if the cells are local. */
        if (ci->nodeID == engine_rank) cell_activate_drift_spart(ci, s);
        if (cj->nodeID == engine_rank) cell_activate_drift_part(cj, s);
        if (cj->nodeID == engine_rank && with_timestep_sync)
          cell_activate_sync_part(cj, s);

        /* Do we need to sort the cells? */
        cell_activate_hydro_sorts(cj, sid, s);
        cell_activate_stars_sorts(ci, sid, s);
      }

      if (cj_active) {

        /* We are going to interact this pair, so store some values. */
        atomic_or(&cj->stars.requires_sorts, 1 << sid);
        atomic_or(&ci->hydro.requires_sorts, 1 << sid);

        ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;
        cj->stars.dx_max_sort_old = cj->stars.dx_max_sort;

        /* Activate the drifts if the cells are local. */
        if (ci->nodeID == engine_rank) cell_activate_drift_part(ci, s);
        if (cj->nodeID == engine_rank) cell_activate_drift_spart(cj, s);
        if (ci->nodeID == engine_rank && with_timestep_sync)
          cell_activate_sync_part(ci, s);

        /* Do we need to sort the cells? */
        cell_activate_hydro_sorts(ci, sid, s);
        cell_activate_stars_sorts(cj, sid, s);
      }
    }
  } /* Otherwise, pair interation */
}

/**
 * @brief Traverse a sub-cell task and activate the black_holes drift tasks that
 * are required by a black_holes task
 *
 * @param ci The first #cell we recurse in.
 * @param cj The second #cell we recurse in.
 * @param s The task #scheduler.
 * @param with_timestep_sync Are we running with time-step synchronization on?
 */
void cell_activate_subcell_black_holes_tasks(struct cell *ci, struct cell *cj,
                                             struct scheduler *s,
                                             const int with_timestep_sync) {
  const struct engine *e = s->space->e;

  /* Store the current dx_max and h_max values. */
  ci->black_holes.dx_max_part_old = ci->black_holes.dx_max_part;
  ci->black_holes.h_max_old = ci->black_holes.h_max;
  ci->hydro.dx_max_part_old = ci->hydro.dx_max_part;
  ci->hydro.h_max_old = ci->hydro.h_max;

  if (cj != NULL) {
    cj->black_holes.dx_max_part_old = cj->black_holes.dx_max_part;
    cj->black_holes.h_max_old = cj->black_holes.h_max;
    cj->hydro.dx_max_part_old = cj->hydro.dx_max_part;
    cj->hydro.h_max_old = cj->hydro.h_max;
  }

  /* Self interaction? */
  if (cj == NULL) {
    /* Do anything? */
    if (!cell_is_active_black_holes(ci, e) || ci->hydro.count == 0 ||
        ci->black_holes.count == 0)
      return;

    /* Recurse? */
    if (cell_can_recurse_in_self_black_holes_task(ci)) {
      /* Loop over all progenies and pairs of progenies */
      for (int j = 0; j < 8; j++) {
        if (ci->progeny[j] != NULL) {
          cell_activate_subcell_black_holes_tasks(ci->progeny[j], NULL, s,
                                                  with_timestep_sync);
          for (int k = j + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              cell_activate_subcell_black_holes_tasks(
                  ci->progeny[j], ci->progeny[k], s, with_timestep_sync);
        }
      }
    } else {
      /* We have reached the bottom of the tree: activate drift */
      cell_activate_drift_bpart(ci, s);
      cell_activate_drift_part(ci, s);
      if (with_timestep_sync) cell_activate_sync_part(ci, s);
    }
  }

  /* Otherwise, pair interation */
  else {

    const int ci_active = cell_is_active_black_holes(ci, e);
    const int cj_active = cell_is_active_black_holes(cj, e);

    /* Should we even bother? */
    if (!ci_active && !cj_active) return;

    /* Get the orientation of the pair. */
    double shift[3];
    const int sid = space_getsid(s->space, &ci, &cj, shift);

    /* recurse? */
    if (cell_can_recurse_in_pair_black_holes_task(ci, cj) &&
        cell_can_recurse_in_pair_black_holes_task(cj, ci)) {
      const struct cell_split_pair *csp = &cell_split_pairs[sid];
      for (int k = 0; k < csp->count; k++) {
        const int pid = csp->pairs[k].pid;
        const int pjd = csp->pairs[k].pjd;
        if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
          cell_activate_subcell_black_holes_tasks(
              ci->progeny[pid], cj->progeny[pjd], s, with_timestep_sync);
      }
    }

    /* Otherwise, activate the drifts. */
    else if (ci_active || cj_active) {

      /* Note we need to drift *both* BH cells to deal with BH<->BH swallows
       * But we only need to drift the gas cell if the *other* cell has an
       * active BH */

      /* Activate the drifts if the cells are local. */
      if (ci->nodeID == engine_rank) cell_activate_drift_bpart(ci, s);
      if (cj->nodeID == engine_rank && ci_active)
        cell_activate_drift_part(cj, s);
      if (cj->nodeID == engine_rank && ci_active && with_timestep_sync)
        cell_activate_sync_part(cj, s);

      /* Activate the drifts if the cells are local. */
      if (ci->nodeID == engine_rank && cj_active)
        cell_activate_drift_part(ci, s);
      if (cj->nodeID == engine_rank) cell_activate_drift_bpart(cj, s);
      if (ci->nodeID == engine_rank && cj_active && with_timestep_sync)
        cell_activate_sync_part(ci, s);
    }
  } /* Otherwise, pair interation */
}

/**
 * @brief Traverse a sub-cell task and activate the sinks drift tasks that
 * are required by a sinks task
 *
 * @param ci The first #cell we recurse in.
 * @param cj The second #cell we recurse in.
 * @param s The task #scheduler.
 * @param with_timestep_sync Are we running with time-step synchronization on?
 */
void cell_activate_subcell_sinks_tasks(struct cell *ci, struct cell *cj,
                                       struct scheduler *s,
                                       const int with_timestep_sync) {
  const struct engine *e = s->space->e;

  /* Store the current dx_max and h_max values. */
  ci->sinks.dx_max_part_old = ci->sinks.dx_max_part;
  ci->sinks.r_cut_max_old = ci->sinks.r_cut_max;
  ci->hydro.dx_max_part_old = ci->hydro.dx_max_part;
  ci->hydro.h_max_old = ci->hydro.h_max;

  if (cj != NULL) {
    cj->sinks.dx_max_part_old = cj->sinks.dx_max_part;
    cj->sinks.r_cut_max_old = cj->sinks.r_cut_max;
    cj->hydro.dx_max_part_old = cj->hydro.dx_max_part;
    cj->hydro.h_max_old = cj->hydro.h_max;
  }

  /* Self interaction? */
  if (cj == NULL) {

    const int ci_active =
        cell_is_active_sinks(ci, e) || cell_is_active_hydro(ci, e);

    /* Do anything? */
    if (!ci_active || (ci->hydro.count == 0 && ci->sinks.count == 0)) return;

    /* Recurse? */
    if (cell_can_recurse_in_self_sinks_task(ci)) {
      /* Loop over all progenies and pairs of progenies */
      for (int j = 0; j < 8; j++) {
        if (ci->progeny[j] != NULL) {
          cell_activate_subcell_sinks_tasks(ci->progeny[j], NULL, s,
                                            with_timestep_sync);
          for (int k = j + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              cell_activate_subcell_sinks_tasks(ci->progeny[j], ci->progeny[k],
                                                s, with_timestep_sync);
        }
      }
    } else {
      /* We have reached the bottom of the tree: activate drift */
      cell_activate_drift_sink(ci, s);
      cell_activate_drift_part(ci, s);
      cell_activate_sink_formation_tasks(ci->top, s);
      if (with_timestep_sync) cell_activate_sync_part(ci, s);
    }
  }

  /* Otherwise, pair interation */
  else {
    /* Get the orientation of the pair. */
    double shift[3];
    const int sid = space_getsid(s->space, &ci, &cj, shift);

    const int ci_active =
        cell_is_active_sinks(ci, e) || cell_is_active_hydro(ci, e);
    const int cj_active =
        cell_is_active_sinks(cj, e) || cell_is_active_hydro(cj, e);

    /* Should we even bother? */
    if (!ci_active && !cj_active) return;

    /* recurse? */
    if (cell_can_recurse_in_pair_sinks_task(ci, cj) &&
        cell_can_recurse_in_pair_sinks_task(cj, ci)) {

      const struct cell_split_pair *csp = &cell_split_pairs[sid];
      for (int k = 0; k < csp->count; k++) {
        const int pid = csp->pairs[k].pid;
        const int pjd = csp->pairs[k].pjd;
        if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
          cell_activate_subcell_sinks_tasks(ci->progeny[pid], cj->progeny[pjd],
                                            s, with_timestep_sync);
      }
    }

    /* Otherwise, activate the sorts and drifts. */
    else {

      /* For the sink mergers */
      if (ci->nodeID == engine_rank) {
        cell_activate_drift_sink(ci, s);
        cell_activate_sink_formation_tasks(ci->top, s);
      }
      if (cj->nodeID == engine_rank) {
        cell_activate_drift_sink(cj, s);
        if (ci->top != cj->top) {
          cell_activate_sink_formation_tasks(cj->top, s);
        }
      }

      if (ci_active) {

        /* We are going to interact this pair, so store some values. */
        atomic_or(&cj->hydro.requires_sorts, 1 << sid);
        cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;

        /* Activate the drifts if the cells are local. */
        if (cj->nodeID == engine_rank) cell_activate_drift_part(cj, s);
        if (cj->nodeID == engine_rank && with_timestep_sync)
          cell_activate_sync_part(cj, s);

        /* Do we need to sort the cells? */
        cell_activate_hydro_sorts(cj, sid, s);
      }

      if (cj_active) {

        /* We are going to interact this pair, so store some values. */
        atomic_or(&ci->hydro.requires_sorts, 1 << sid);
        ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;

        /* Activate the drifts if the cells are local. */
        if (ci->nodeID == engine_rank) cell_activate_drift_part(ci, s);
        if (ci->nodeID == engine_rank && with_timestep_sync)
          cell_activate_sync_part(ci, s);

        /* Do we need to sort the cells? */
        cell_activate_hydro_sorts(ci, sid, s);
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
 *
 * @return 1 if the cell and its progeny have fully been treated.
 */
int cell_activate_subcell_grav_tasks(struct cell *restrict ci,
                                     struct cell *restrict cj,
                                     struct scheduler *s) {

  /* Some constants */
  const struct space *sp = s->space;
  const struct engine *e = sp->e;
  const int nodeID = e->nodeID;

  /* Self interaction? */
  if (cj == NULL) {

    /* Do anything? */
    if (ci->grav.count == 0 || !cell_is_active_gravity(ci, e)) return 1;

    /* Has it already been processed? */
    if (cell_get_flag(ci, cell_flag_unskip_self_grav_processed)) return 1;

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

    /* Flag the cell has having been treated */
    cell_set_flag(ci, cell_flag_unskip_self_grav_processed);

    /* And also return that information */
    return 1;
  }

  /* Pair interaction */
  else {

    /* Has it already been processed? */
    if (cell_get_flag(ci, cell_flag_unskip_pair_grav_processed) &&
        cell_get_flag(cj, cell_flag_unskip_pair_grav_processed))
      return 1;

    /* Anything to do here?
     * Note that Here we return 0 as another pair direction for either of
     * the two cells could pass the tests. */
    const int do_ci = cell_is_active_gravity(ci, e) && ci->nodeID == nodeID;
    const int do_cj = cell_is_active_gravity(cj, e) && cj->nodeID == nodeID;
    if (!do_ci && !do_cj) return 0;
    if (ci->grav.count == 0 || cj->grav.count == 0) return 0;

    /* Atomically drift the multipole in ci */
    lock_lock(&ci->grav.mlock);
    if (ci->grav.ti_old_multipole < e->ti_current) cell_drift_multipole(ci, e);
    if (lock_unlock(&ci->grav.mlock) != 0) error("Impossible to unlock m-pole");

    /* Atomically drift the multipole in cj */
    lock_lock(&cj->grav.mlock);
    if (cj->grav.ti_old_multipole < e->ti_current) cell_drift_multipole(cj, e);
    if (lock_unlock(&cj->grav.mlock) != 0) error("Impossible to unlock m-pole");

    /* Can we use multipoles ? */
    if (cell_can_use_pair_mm(ci, cj, e, sp, /*use_rebuild_data=*/0,
                             /*is_tree_walk=*/1)) {

      /* Ok, no need to drift anything */
      return 0;
    }

    /* Otherwise, if we are at the bottom, activate the gpart drifts. */
    else if (!ci->split && !cj->split) {

      /* Activate the drifts if the cells are local. */
      if (cell_is_active_gravity(ci, e) || cell_is_active_gravity(cj, e)) {

        cell_set_flag(ci, cell_flag_unskip_pair_grav_processed);
        cell_set_flag(cj, cell_flag_unskip_pair_grav_processed);

        if (ci->nodeID == nodeID) cell_activate_drift_gpart(ci, s);
        if (cj->nodeID == nodeID) cell_activate_drift_gpart(cj, s);
      }

      /* And return that information */
      return 1;

    }

    /* Ok, we can still recurse at least on one side. */
    else {

      /* Recover the multipole information */
      const struct gravity_tensors *const multi_i = ci->grav.multipole;
      const struct gravity_tensors *const multi_j = cj->grav.multipole;
      const double ri_max = multi_i->r_max;
      const double rj_max = multi_j->r_max;

      int ci_number_children = 0, cj_number_children = 0;
      int progenies_all_processed = 0;

      /* Let's open up the largest of the two cells,
       * provided it is split into smaller cells. */

      if (ri_max > rj_max) {

        if (ci->split) {

          /* Loop over ci's children, activate what is needed and
           * collect the number of cells that have been fully processed */
          for (int k = 0; k < 8; k++) {
            if (ci->progeny[k] != NULL) {
              ci_number_children++;
              progenies_all_processed +=
                  cell_activate_subcell_grav_tasks(ci->progeny[k], cj, s);
            }
          }

          const int cell_done = (progenies_all_processed == ci_number_children);

          /* Flag the cells as being fully processed. */
          if (cell_done) {
            cell_set_flag(ci, cell_flag_unskip_pair_grav_processed);
          }

          return cell_done;

        } else if (cj->split) {

          /* Loop over cj's children, activate what is needed and
           * collect the number of cells that have been fully processed */
          for (int k = 0; k < 8; k++) {
            if (cj->progeny[k] != NULL) {
              cj_number_children++;
              progenies_all_processed +=
                  cell_activate_subcell_grav_tasks(ci, cj->progeny[k], s);
            }
          }

          const int cell_done = (progenies_all_processed == cj_number_children);

          /* Flag the cells as being fully processed. */
          if (cell_done) {
            cell_set_flag(cj, cell_flag_unskip_pair_grav_processed);
          }

          return cell_done;

        } else {

#ifdef SWIFT_DEBUG_CHECKS
          error("Fundamental error in the logic");
#endif
        }

      } else if (rj_max >= ri_max) {

        if (cj->split) {

          /* Loop over cj's children, activate what is needed and
           * collect the number of cells that have been fully processed */
          for (int k = 0; k < 8; k++) {
            if (cj->progeny[k] != NULL) {
              cj_number_children++;
              progenies_all_processed +=
                  cell_activate_subcell_grav_tasks(ci, cj->progeny[k], s);
            }
          }

          const int cell_done = (progenies_all_processed == cj_number_children);

          /* Flag the cells as being fully processed. */
          if (cell_done) {
            cell_set_flag(cj, cell_flag_unskip_pair_grav_processed);
          }

          return cell_done;

        } else if (ci->split) {

          /* Loop over ci's children, activate what is needed and
           * collect the number of cells that have been fully processed */
          for (int k = 0; k < 8; k++) {
            if (ci->progeny[k] != NULL) {
              ci_number_children++;
              progenies_all_processed +=
                  cell_activate_subcell_grav_tasks(ci->progeny[k], cj, s);
            }
          }

          const int cell_done = (progenies_all_processed == ci_number_children);

          /* Flag the cells as being done. */
          if (cell_done) {
            cell_set_flag(ci, cell_flag_unskip_pair_grav_processed);
          }

          return cell_done;

        } else {
#ifdef SWIFT_DEBUG_CHECKS
          error("Fundamental error in the logic");
#endif
        }
      }
    }
  }
#ifdef SWIFT_DEBUG_CHECKS
  error("Fundamental error in the logic");
#endif
  return -1;
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
 * @brief Traverse a sub-cell task and activate the sort tasks that are
 * required by a RT task
 *
 * @param ci The first #cell we recurse in.
 * @param cj The second #cell we recurse in.
 * @param s The task #scheduler.
 * @param sub_cycle Are we in a subcycle or not?
 */
void cell_activate_subcell_rt_tasks(struct cell *ci, struct cell *cj,
                                    struct scheduler *s, const int sub_cycle) {

  /* Only do this during real time steps, not during subcycling. */
  if (sub_cycle) return;
  const struct engine *e = s->space->e;

  /* Store the current dx_max and h_max values. */
  ci->hydro.dx_max_part_old = ci->hydro.dx_max_part;
  ci->hydro.h_max_old = ci->hydro.h_max;

  if (cj != NULL) {
    cj->hydro.dx_max_part_old = cj->hydro.dx_max_part;
    cj->hydro.h_max_old = cj->hydro.h_max;
  }

  const int ci_active = cell_is_rt_active(ci, e);
  const int cj_active = ((cj != NULL) && cell_is_rt_active(cj, e));

  /* Self interaction? */
  if (cj == NULL) {
    /* Do anything? */
    if (ci->hydro.count == 0 || !ci_active) return;

    /* Recurse? */
    if (cell_can_recurse_in_self_hydro_task(ci)) {
      /* Loop over all progenies and pairs of progenies */
      for (int j = 0; j < 8; j++) {
        if (ci->progeny[j] != NULL) {
          cell_activate_subcell_rt_tasks(ci->progeny[j], NULL, s, sub_cycle);
          for (int k = j + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              cell_activate_subcell_rt_tasks(ci->progeny[j], ci->progeny[k], s,
                                             sub_cycle);
        }
      }
    }
  }

  /* Otherwise, pair interation */
  else {

    /* Should we even bother? */
    if (!ci_active && !cj_active) return;
    if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

    /* Get the orientation of the pair. */
    double shift[3];
    const int sid = space_getsid(s->space, &ci, &cj, shift);

    /* recurse? */
    if (cell_can_recurse_in_pair_hydro_task(ci) &&
        cell_can_recurse_in_pair_hydro_task(cj)) {
      const struct cell_split_pair *csp = &cell_split_pairs[sid];
      for (int k = 0; k < csp->count; k++) {
        const int pid = csp->pairs[k].pid;
        const int pjd = csp->pairs[k].pjd;
        if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
          cell_activate_subcell_rt_tasks(ci->progeny[pid], cj->progeny[pjd], s,
                                         sub_cycle);
      }
    }

    /* Otherwise, activate the sorts and drifts. */
    else if (ci_active || cj_active) {

      /* We are going to interact this pair, so store some values. */
      atomic_or(&ci->hydro.requires_sorts, 1 << sid);
      atomic_or(&cj->hydro.requires_sorts, 1 << sid);
      ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;
      cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;

      /* Do we need to sort the cells? */
      cell_activate_rt_sorts(ci, sid, s);
      cell_activate_rt_sorts(cj, sid, s);
    }
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
  const int with_feedback = e->policy & engine_policy_feedback;
  const int with_timestep_limiter =
      (e->policy & engine_policy_timestep_limiter);
  const int with_sinks = e->policy & engine_policy_sinks;

#ifdef WITH_MPI
  const int with_star_formation = e->policy & engine_policy_star_formation;
  if (with_sinks) error("Cannot use sink tasks and MPI");
#endif
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
        if (ci_nodeID == nodeID && with_timestep_limiter)
          cell_activate_limiter(ci, s);
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
        if (ci_nodeID == nodeID && with_timestep_limiter)
          cell_activate_limiter(ci, s);
        if (cj_nodeID == nodeID && with_timestep_limiter)
          cell_activate_limiter(cj, s);

        /* Check the sorts and activate them if needed. */
        cell_activate_hydro_sorts(ci, t->flags, s);
        cell_activate_hydro_sorts(cj, t->flags, s);
      }

      /* Store current values of dx_max and h_max. */
      else if (t->type == task_type_sub_self) {
        cell_activate_subcell_hydro_tasks(ci, NULL, s, with_timestep_limiter);
      }

      /* Store current values of dx_max and h_max. */
      else if (t->type == task_type_sub_pair) {
        cell_activate_subcell_hydro_tasks(ci, cj, s, with_timestep_limiter);
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
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_xv);
          if (ci_active) {
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_rho);

#ifdef EXTRA_HYDRO_LOOP
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_gradient);
#endif
          }
        }
        /* If the local cell is inactive and the remote cell is active, we
         * still need to receive stuff to be able to do the force interaction
         * on this node as well. */
        else if (ci_active) {
#ifdef MPI_SYMMETRIC_FORCE_INTERACTION
          /* NOTE: (yuyttenh, 09/2022) Since the particle communications send
           * over whole particles currently, just activating the gradient
           * send/recieve should be enough for now. The remote active
           * particles are only needed for the sorts and the flux exchange on
           * the node of the inactive cell, so sending over the xv and
           * gradient suffices. If at any point the commutications change, we
           * should probably also send over the rho separately. */
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_xv);
#ifndef EXTRA_HYDRO_LOOP
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_rho);
#else
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_gradient);
#endif
#endif
        }

        /* If the foreign cell is active, we want its particles for the limiter
         */
        if (ci_active && with_timestep_limiter) {
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_limiter);
          scheduler_activate_unpack(s, ci->mpi.unpack, task_subtype_limiter);
        }

        /* Is the foreign cell active and will need stuff from us? */
        if (ci_active) {

          scheduler_activate_send(s, cj->mpi.send, task_subtype_xv, ci_nodeID);

          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(cj, s);
          if (with_timestep_limiter) cell_activate_limiter(cj, s);

          /* If the local cell is also active, more stuff will be needed. */
          if (cj_active) {
            scheduler_activate_send(s, cj->mpi.send, task_subtype_rho,
                                    ci_nodeID);

#ifdef EXTRA_HYDRO_LOOP
            scheduler_activate_send(s, cj->mpi.send, task_subtype_gradient,
                                    ci_nodeID);
#endif
          }
        }
        /* If the foreign cell is inactive, but the local cell is active,
         * we still need to send stuff to be able to do the force interaction
         * on both nodes */
        else if (cj_active) {
#ifdef MPI_SYMMETRIC_FORCE_INTERACTION
          /* See NOTE on line 1542 */
          scheduler_activate_send(s, cj->mpi.send, task_subtype_xv, ci_nodeID);
          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(cj, s);
#ifndef EXTRA_HYDRO_LOOP
          scheduler_activate_send(s, cj->mpi.send, task_subtype_rho, ci_nodeID);
#else
          scheduler_activate_send(s, cj->mpi.send, task_subtype_gradient,
                                  ci_nodeID);
#endif
#endif
        }

        /* If the local cell is active, send its particles for the limiting. */
        if (cj_active && with_timestep_limiter) {
          scheduler_activate_send(s, cj->mpi.send, task_subtype_limiter,
                                  ci_nodeID);
          scheduler_activate_pack(s, cj->mpi.pack, task_subtype_limiter,
                                  ci_nodeID);
        }

        /* Propagating new star counts? */
        if (with_star_formation && with_feedback) {
          if (ci_active && ci->hydro.count > 0) {
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_sf_counts);
          }
          if (cj_active && cj->hydro.count > 0) {
            scheduler_activate_send(s, cj->mpi.send, task_subtype_sf_counts,
                                    ci_nodeID);
          }
        }

      } else if (cj_nodeID != nodeID) {
        /* If the local cell is active, receive data from the foreign cell. */
        if (ci_active) {
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_xv);
          if (cj_active) {
            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_rho);

#ifdef EXTRA_HYDRO_LOOP
            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_gradient);
#endif
          }
        }
        /* If the local cell is inactive and the remote cell is active, we
         * still need to receive stuff to be able to do the force interaction
         * on this node as well. */
        else if (cj_active) {
#ifdef MPI_SYMMETRIC_FORCE_INTERACTION
          /* See NOTE on line 1542. */
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_xv);
#ifndef EXTRA_HYDRO_LOOP
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_rho);
#else
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_gradient);
#endif
#endif
        }

        /* If the foreign cell is active, we want its particles for the limiter
         */
        if (cj_active && with_timestep_limiter) {
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_limiter);
          scheduler_activate_unpack(s, cj->mpi.unpack, task_subtype_limiter);
        }

        /* Is the foreign cell active and will need stuff from us? */
        if (cj_active) {

          scheduler_activate_send(s, ci->mpi.send, task_subtype_xv, cj_nodeID);

          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(ci, s);
          if (with_timestep_limiter) cell_activate_limiter(ci, s);

          /* If the local cell is also active, more stuff will be needed. */
          if (ci_active) {

            scheduler_activate_send(s, ci->mpi.send, task_subtype_rho,
                                    cj_nodeID);

#ifdef EXTRA_HYDRO_LOOP
            scheduler_activate_send(s, ci->mpi.send, task_subtype_gradient,
                                    cj_nodeID);
#endif
          }
        }
        /* If the foreign cell is inactive, but the local cell is active,
         * we still need to send stuff to be able to do the force interaction
         * on both nodes */
        else if (ci_active) {
#ifdef MPI_SYMMETRIC_FORCE_INTERACTION
          /* See NOTE on line 1542. */
          scheduler_activate_send(s, ci->mpi.send, task_subtype_xv, cj_nodeID);
          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(ci, s);
#ifndef EXTRA_HYDRO_LOOP
          scheduler_activate_send(s, ci->mpi.send, task_subtype_rho, cj_nodeID);
#else
          scheduler_activate_send(s, ci->mpi.send, task_subtype_gradient,
                                  cj_nodeID);
#endif
#endif
        }

        /* If the local cell is active, send its particles for the limiting. */
        if (ci_active && with_timestep_limiter) {
          scheduler_activate_send(s, ci->mpi.send, task_subtype_limiter,
                                  cj_nodeID);
          scheduler_activate_pack(s, ci->mpi.pack, task_subtype_limiter,
                                  cj_nodeID);
        }

        /* Propagating new star counts? */
        if (with_star_formation && with_feedback) {
          if (cj_active && cj->hydro.count > 0) {
            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_sf_counts);
          }
          if (ci_active && ci->hydro.count > 0) {
            scheduler_activate_send(s, ci->mpi.send, task_subtype_sf_counts,
                                    cj_nodeID);
          }
        }
      }
#endif
    }
  }

  /* Unskip all the other task types. */
  int c_active = cell_is_active_hydro(c, e);
  if (c->nodeID == nodeID && c_active) {
    for (struct link *l = c->hydro.gradient; l != NULL; l = l->next) {
      scheduler_activate(s, l->t);
    }
    for (struct link *l = c->hydro.force; l != NULL; l = l->next) {
      scheduler_activate(s, l->t);
    }

    for (struct link *l = c->hydro.limiter; l != NULL; l = l->next)
      scheduler_activate(s, l->t);

    if (c->hydro.extra_ghost != NULL)
      scheduler_activate(s, c->hydro.extra_ghost);
    if (c->hydro.ghost_in != NULL) cell_activate_hydro_ghosts(c, s, e);
    if (c->kick1 != NULL) scheduler_activate(s, c->kick1);
    if (c->kick2 != NULL) scheduler_activate(s, c->kick2);
    if (c->timestep != NULL) scheduler_activate(s, c->timestep);
    if (c->top->timestep_collect != NULL)
      scheduler_activate(s, c->top->timestep_collect);
    if (c->hydro.end_force != NULL) scheduler_activate(s, c->hydro.end_force);
    if (c->hydro.cooling_in != NULL) cell_activate_cooling(c, s, e);
#ifdef WITH_CSDS
    if (c->csds != NULL) scheduler_activate(s, c->csds);
#endif

    if (c->top->hydro.star_formation != NULL) {
      cell_activate_star_formation_tasks(c->top, s, with_feedback);
      cell_activate_super_spart_drifts(c->top, s);
    }
    if (with_sinks && c->top->sinks.star_formation_sink != NULL) {
      cell_activate_star_formation_sink_tasks(c->top, s, with_feedback);
    }
  }
  /* Additionally unskip force interactions between inactive local cell and
   * active remote cell. (The cell unskip will only be called for active cells,
   * so, we have to do this now, from the active remote cell). */
  else if (c->nodeID != nodeID && c_active) {
#if defined(MPI_SYMMETRIC_FORCE_INTERACTION) && defined(WITH_MPI)
    for (struct link *l = c->hydro.force; l != NULL; l = l->next) {
      struct task *t = l->t;
      if (t->type != task_type_pair && t->type != task_type_sub_pair) continue;

      struct cell *ci = l->t->ci;
      struct cell *cj = l->t->cj;

      const int ci_active = cell_is_active_hydro(ci, e);
      const int cj_active = cell_is_active_hydro(cj, e);
      const int ci_nodeID = ci->nodeID;
      const int cj_nodeID = cj->nodeID;
      if ((!ci_active && ci_nodeID == nodeID && cj_active &&
           cj_nodeID != nodeID) ||
          (!cj_active && cj_nodeID == nodeID && ci_active &&
           ci_nodeID != nodeID)) {
        scheduler_activate(s, l->t);

        if (t->type == task_type_pair) {
          /* Store some values. */
          atomic_or(&ci->hydro.requires_sorts, 1 << t->flags);
          atomic_or(&cj->hydro.requires_sorts, 1 << t->flags);
          ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;
          cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;

          /* Activate the drift tasks. */
          if (ci_nodeID == nodeID) cell_activate_drift_part(ci, s);
          if (cj_nodeID == nodeID) cell_activate_drift_part(cj, s);

          /* Activate the limiter tasks. */
          if (ci_nodeID == nodeID && with_timestep_limiter)
            cell_activate_limiter(ci, s);
          if (cj_nodeID == nodeID && with_timestep_limiter)
            cell_activate_limiter(cj, s);

          /* Check the sorts and activate them if needed. */
          cell_activate_hydro_sorts(ci, t->flags, s);
          cell_activate_hydro_sorts(cj, t->flags, s);
        }

        /* Store current values of dx_max and h_max. */
        else if (t->type == task_type_sub_pair) {
          cell_activate_subcell_hydro_tasks(ci, cj, s, with_timestep_limiter);
        }
      }
    }
#endif
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
        if (cj_active)
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_gpart);

        /* Is the foreign cell active and will need stuff from us? */
        if (ci_active) {

          scheduler_activate_send(s, cj->mpi.send, task_subtype_gpart,
                                  ci_nodeID);

          /* Drift the cell which will be sent at the level at which it is
             sent, i.e. drift the cell specified in the send task (l->t)
             itself. */
          cell_activate_drift_gpart(cj, s);
        }

      } else if (cj_nodeID != nodeID) {
        /* If the local cell is active, receive data from the foreign cell. */
        if (ci_active)
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_gpart);

        /* Is the foreign cell active and will need stuff from us? */
        if (cj_active) {

          scheduler_activate_send(s, ci->mpi.send, task_subtype_gpart,
                                  cj_nodeID);

          /* Drift the cell which will be sent at the level at which it is
             sent, i.e. drift the cell specified in the send task (l->t)
             itself. */
          cell_activate_drift_gpart(ci, s);
        }
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
    if (c->top->timestep_collect != NULL)
      scheduler_activate(s, c->top->timestep_collect);
    if (c->grav.down != NULL) scheduler_activate(s, c->grav.down);
    if (c->grav.down_in != NULL) scheduler_activate(s, c->grav.down_in);
    if (c->grav.long_range != NULL) scheduler_activate(s, c->grav.long_range);
    if (c->grav.end_force != NULL) scheduler_activate(s, c->grav.end_force);
    if (c->grav.neutrino_weight != NULL)
      scheduler_activate(s, c->grav.neutrino_weight);
#ifdef WITH_CSDS
    if (c->csds != NULL) scheduler_activate(s, c->csds);
#endif
  }

  return rebuild;
}

/**
 * @brief Un-skips all the stars tasks associated with a given cell and checks
 * if the space needs to be rebuilt.
 *
 * @param c the #cell.
 * @param s the #scheduler.
 * @param with_star_formation Are we running with star formation switched on?
 * @param with_star_formation_sink Are we running with star formation based on
 * sink switched on?
 *
 * @return 1 If the space needs rebuilding. 0 otherwise.
 */
int cell_unskip_stars_tasks(struct cell *c, struct scheduler *s,
                            const int with_star_formation,
                            const int with_star_formation_sink) {

  struct engine *e = s->space->e;
  const int with_timestep_sync = (e->policy & engine_policy_timestep_sync);
  const int nodeID = e->nodeID;
  int rebuild = 0;

  if (c->stars.drift != NULL) {
    if (cell_need_activating_stars(c, e, with_star_formation,
                                   with_star_formation_sink)) {

      cell_activate_drift_spart(c, s);
    }
  }

  /* Un-skip the density tasks involved with this cell. */
  for (struct link *l = c->stars.density; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    const int ci_active = cell_need_activating_stars(ci, e, with_star_formation,
                                                     with_star_formation_sink);

    const int cj_active =
        (cj != NULL) && cell_need_activating_stars(cj, e, with_star_formation,
                                                   with_star_formation_sink);

    /* Activate the drifts */
    if (t->type == task_type_self && ci_active) {
      cell_activate_drift_spart(ci, s);
      cell_activate_drift_part(ci, s);
      if (with_timestep_sync) cell_activate_sync_part(ci, s);
    }

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active || cj_active) &&
        (ci_nodeID == nodeID || cj_nodeID == nodeID)) {
      scheduler_activate(s, t);

      if (t->type == task_type_pair) {
        /* Activate stars_in for each cell that is part of
         * a pair as to not miss any dependencies */
        if (ci_nodeID == nodeID)
          scheduler_activate(s, ci->hydro.super->stars.stars_in);
        if (cj_nodeID == nodeID)
          scheduler_activate(s, cj->hydro.super->stars.stars_in);

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
          if (cj_nodeID == nodeID && with_timestep_sync)
            cell_activate_sync_part(cj, s);

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
          if (ci_nodeID == nodeID && with_timestep_sync)
            cell_activate_sync_part(ci, s);

          /* Check the sorts and activate them if needed. */
          cell_activate_hydro_sorts(ci, t->flags, s);
          cell_activate_stars_sorts(cj, t->flags, s);
        }
      }

      else if (t->type == task_type_sub_self) {
        cell_activate_subcell_stars_tasks(ci, NULL, s, with_star_formation,
                                          with_star_formation_sink,
                                          with_timestep_sync);
      }

      else if (t->type == task_type_sub_pair) {
        cell_activate_subcell_stars_tasks(ci, cj, s, with_star_formation,
                                          with_star_formation_sink,
                                          with_timestep_sync);

        /* Activate stars_in for each cell that is part of
         * a sub_pair task as to not miss any dependencies */
        if (ci_nodeID == nodeID)
          scheduler_activate(s, ci->hydro.super->stars.stars_in);
        if (cj_nodeID == nodeID)
          scheduler_activate(s, cj->hydro.super->stars.stars_in);
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
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_xv);
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_rho);
#ifdef EXTRA_STAR_LOOPS
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_part_prep1);
#endif
          /* If the local cell is active, more stuff will be needed. */
          scheduler_activate_send(s, cj->mpi.send, task_subtype_spart_density,
                                  ci_nodeID);
#ifdef EXTRA_STAR_LOOPS
          scheduler_activate_send(s, cj->mpi.send, task_subtype_spart_prep2,
                                  ci_nodeID);
#endif
          cell_activate_drift_spart(cj, s);
        }

        if (ci_active) {
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_spart_density);
#ifdef EXTRA_STAR_LOOPS
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_spart_prep2);
#endif

          /* Is the foreign cell active and will need stuff from us? */
          scheduler_activate_send(s, cj->mpi.send, task_subtype_xv, ci_nodeID);
          scheduler_activate_send(s, cj->mpi.send, task_subtype_rho, ci_nodeID);
#ifdef EXTRA_STAR_LOOPS
          scheduler_activate_send(s, cj->mpi.send, task_subtype_part_prep1,
                                  ci_nodeID);
#endif
          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(cj, s);
        }

      } else if (cj_nodeID != nodeID) {
        /* If the local cell is active, receive data from the foreign cell. */
        if (ci_active) {
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_xv);
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_rho);
#ifdef EXTRA_STAR_LOOPS
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_part_prep1);
#endif
          /* If the local cell is active, more stuff will be needed. */
          scheduler_activate_send(s, ci->mpi.send, task_subtype_spart_density,
                                  cj_nodeID);
#ifdef EXTRA_STAR_LOOPS
          scheduler_activate_send(s, ci->mpi.send, task_subtype_spart_prep2,
                                  cj_nodeID);
#endif
          cell_activate_drift_spart(ci, s);
        }

        if (cj_active) {
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_spart_density);
#ifdef EXTRA_STAR_LOOPS
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_spart_prep2);
#endif

          /* Is the foreign cell active and will need stuff from us? */
          scheduler_activate_send(s, ci->mpi.send, task_subtype_xv, cj_nodeID);
          scheduler_activate_send(s, ci->mpi.send, task_subtype_rho, cj_nodeID);
#ifdef EXTRA_STAR_LOOPS
          scheduler_activate_send(s, ci->mpi.send, task_subtype_part_prep1,
                                  cj_nodeID);
#endif

          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(ci, s);
        }
      }
#endif
    }
  }

  for (struct link *l = c->stars.prepare1; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    const int ci_active = cell_need_activating_stars(ci, e, with_star_formation,
                                                     with_star_formation_sink);

    const int cj_active =
        (cj != NULL) && cell_need_activating_stars(cj, e, with_star_formation,
                                                   with_star_formation_sink);

#ifdef SWIFT_DEBUG_CHECKS
    if (with_star_formation_sink) {
      error("TODO");
    }
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

        /* If there are active sparts in cj, activate hydro ghost in ci */
        if (cj_active) {
          scheduler_activate(s, ci->hydro.super->hydro.prep1_ghost);
        }
        /* If there are active sparts in ci, activate hydro ghost in cj */
        else {
          scheduler_activate(s, cj->hydro.super->hydro.prep1_ghost);
        }
      }
      /* Cells ci and cj are from different MPI domains */
      else if ((ci_nodeID == nodeID && cj_nodeID != nodeID) && (cj_active)) {
        /* In task prepare1, we update gas so sparts must be on foreign node */
        scheduler_activate(s, t);
        /* If there are active sparts in cj, activate hydro ghost in ci */
        scheduler_activate(s, ci->hydro.super->hydro.prep1_ghost);
      } else if ((ci_nodeID != nodeID && cj_nodeID == nodeID) && (ci_active)) {
        /* In task prepare1, we update gas so sparts must be on foreign node */
        scheduler_activate(s, t);
        /* If there are active sparts in ci, activate hydro ghost in cj */
        scheduler_activate(s, cj->hydro.super->hydro.prep1_ghost);
      }
    }
  }

  for (struct link *l = c->stars.prepare2; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

#ifdef SWIFT_DEBUG_CHECKS
    if (with_star_formation_sink) {
      error("TODO");
    }
#endif

    const int ci_active = cell_need_activating_stars(ci, e, with_star_formation,
                                                     with_star_formation_sink);

    const int cj_active =
        (cj != NULL) && cell_need_activating_stars(cj, e, with_star_formation,
                                                   with_star_formation_sink);

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
      }
      /* Cells ci and cj are from different MPI domains */
      else if ((ci_nodeID == nodeID && cj_nodeID != nodeID) && (ci_active)) {
        /* In task prep2, we update stars so sparts must be on the local node */
        scheduler_activate(s, t);
      } else if ((ci_nodeID != nodeID && cj_nodeID == nodeID) && (cj_active)) {
        /* In task prep2, we update stars so sparts must be on the local node */
        scheduler_activate(s, t);
      }
    }
  }

  /* Un-skip the feedback tasks involved with this cell. */
  for (struct link *l = c->stars.feedback; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    const int ci_active = cell_need_activating_stars(ci, e, with_star_formation,
                                                     with_star_formation_sink);

    const int cj_active =
        (cj != NULL) && cell_need_activating_stars(cj, e, with_star_formation,
                                                   with_star_formation_sink);

    if (t->type == task_type_self && ci_active) {
      scheduler_activate(s, t);
    }

    else if (t->type == task_type_sub_self && ci_active) {
      scheduler_activate(s, t);
    }

    else if (t->type == task_type_pair || t->type == task_type_sub_pair) {

      if (ci_active || cj_active) {
        /* Activate stars_out for each cell that is part of
         * a pair/sub_pair task as to not miss any dependencies */
        if (ci_nodeID == nodeID)
          scheduler_activate(s, ci->hydro.super->stars.stars_out);
        if (cj_nodeID == nodeID)
          scheduler_activate(s, cj->hydro.super->stars.stars_out);
      }

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
  if (c->nodeID == nodeID) {
    if (cell_need_activating_stars(c, e, with_star_formation,
                                   with_star_formation_sink)) {

      if (c->stars.density_ghost != NULL)
        scheduler_activate(s, c->stars.density_ghost);
      if (c->stars.prep1_ghost != NULL)
        scheduler_activate(s, c->stars.prep1_ghost);
      if (c->hydro.prep1_ghost != NULL)
        scheduler_activate(s, c->hydro.prep1_ghost);
      if (c->stars.prep2_ghost != NULL)
        scheduler_activate(s, c->stars.prep2_ghost);
      /* If we don't have pair tasks, then the stars_in and stars_out still
       * need reactivation. */
      if (c->stars.stars_in != NULL) scheduler_activate(s, c->stars.stars_in);
      if (c->stars.stars_out != NULL) scheduler_activate(s, c->stars.stars_out);
      if (c->kick1 != NULL) scheduler_activate(s, c->kick1);
      if (c->kick2 != NULL) scheduler_activate(s, c->kick2);
      if (c->timestep != NULL) scheduler_activate(s, c->timestep);
      if (c->top->timestep_collect != NULL)
        scheduler_activate(s, c->top->timestep_collect);
#ifdef WITH_CSDS
      if (c->csds != NULL) scheduler_activate(s, c->csds);
#endif
    }
  }

  return rebuild;
}

/**
 * @brief Un-skips all the black hole tasks associated with a given cell and
 * checks if the space needs to be rebuilt.
 *
 * @param c the #cell.
 * @param s the #scheduler.
 *
 * @return 1 If the space needs rebuilding. 0 otherwise.
 */
int cell_unskip_black_holes_tasks(struct cell *c, struct scheduler *s) {

  struct engine *e = s->space->e;
  const int with_timestep_sync = (e->policy & engine_policy_timestep_sync);
  const int nodeID = e->nodeID;
  int rebuild = 0;

  if (c->black_holes.drift != NULL && cell_is_active_black_holes(c, e)) {
    cell_activate_drift_bpart(c, s);
  }

  /* Un-skip the density tasks involved with this cell. */
  for (struct link *l = c->black_holes.density; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const int ci_active =
        ci->black_holes.count > 0 && cell_is_active_black_holes(ci, e);
    const int cj_active = (cj != NULL) ? (cj->black_holes.count > 0 &&
                                          cell_is_active_black_holes(cj, e))
                                       : 0;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active || cj_active) &&
        (ci_nodeID == nodeID || cj_nodeID == nodeID)) {

      scheduler_activate(s, t);

      /* Activate the drifts & sync */
      if (t->type == task_type_self) {
        cell_activate_drift_part(ci, s);
        cell_activate_drift_bpart(ci, s);
        if (with_timestep_sync) cell_activate_sync_part(ci, s);
      }

      /* Activate the drifts */
      else if (t->type == task_type_pair) {

        /* Activate the drift & sync tasks.
         * Note we need to drift *both* BH cells to deal with BH<->BH swallows
         * But we only need to drift the gas cell if the *other* cell has an
         * active BH */
        if (ci_nodeID == nodeID) cell_activate_drift_bpart(ci, s);
        if (ci_nodeID == nodeID && cj_active) cell_activate_drift_part(ci, s);

        if (cj_nodeID == nodeID && ci_active) cell_activate_drift_part(cj, s);
        if (cj_nodeID == nodeID) cell_activate_drift_bpart(cj, s);

        if (ci_nodeID == nodeID && cj_active && with_timestep_sync)
          cell_activate_sync_part(ci, s);
        if (cj_nodeID == nodeID && ci_active && with_timestep_sync)
          cell_activate_sync_part(cj, s);
      }

      /* Store current values of dx_max and h_max. */
      else if (t->type == task_type_sub_self) {
        cell_activate_subcell_black_holes_tasks(ci, NULL, s,
                                                with_timestep_sync);
      }

      /* Store current values of dx_max and h_max. */
      else if (t->type == task_type_sub_pair) {
        cell_activate_subcell_black_holes_tasks(ci, cj, s, with_timestep_sync);
      }

      if (t->type == task_type_pair || t->type == task_type_sub_pair) {
        /* Activate bh_in for each cell that is part of
         * a pair task as to not miss any dependencies */
        if (ci_nodeID == nodeID)
          scheduler_activate(s, ci->hydro.super->black_holes.black_holes_in);
        if (cj_nodeID == nodeID)
          scheduler_activate(s, cj->hydro.super->black_holes.black_holes_in);
      }
    }

    /* Only interested in pair interactions as of here. */
    if (t->type == task_type_pair || t->type == task_type_sub_pair) {

      /* Check whether there was too much particle motion, i.e. the
         cell neighbour conditions were violated. */
      if (cell_need_rebuild_for_black_holes_pair(ci, cj)) rebuild = 1;
      if (cell_need_rebuild_for_black_holes_pair(cj, ci)) rebuild = 1;

      if (ci->hydro.super->black_holes.count > 0 && ci_active)
        scheduler_activate(s, ci->hydro.super->black_holes.swallow_ghost_1);
      if (cj->hydro.super->black_holes.count > 0 && cj_active)
        scheduler_activate(s, cj->hydro.super->black_holes.swallow_ghost_1);

#ifdef WITH_MPI
      /* Activate the send/recv tasks. */
      if (ci_nodeID != nodeID) {

        if (cj_active) {

          /* Receive the foreign parts to compute BH accretion rates and do the
           * swallowing */
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_rho);
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_part_swallow);
          /* scheduler_activate_recv(s, ci->mpi.recv,
           * task_subtype_bpart_merger); */

          /* Send the local BHs to tag the particles to swallow and to do
           * feedback */
          scheduler_activate_send(s, cj->mpi.send, task_subtype_bpart_rho,
                                  ci_nodeID);
          scheduler_activate_send(s, cj->mpi.send, task_subtype_bpart_feedback,
                                  ci_nodeID);
          /* scheduler_activate_send(s, cj->mpi.send, task_subtype_bpart_merger,
           */
          /*                         cj_nodeID); */

          /* Drift before you send */
          cell_activate_drift_bpart(cj, s);
        }

        if (ci_active) {

          /* Receive the foreign BHs to tag particles to swallow and for
           * feedback */
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_bpart_rho);
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_bpart_feedback);

          /* Send the local part information */
          scheduler_activate_send(s, cj->mpi.send, task_subtype_rho, ci_nodeID);
          scheduler_activate_send(s, cj->mpi.send, task_subtype_part_swallow,
                                  ci_nodeID);

          /* if (cj->black_holes.count > 0)  */
          /* scheduler_activate_send(s, cj->mpi.send, task_subtype_bpart_merger,
           */
          /* 			  ci_nodeID); */

          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(cj, s);
          /* if (cj->black_holes.count > 0) cell_activate_drift_bpart(cj, s); */
        }

      } else if (cj_nodeID != nodeID) {

        if (ci_active) {

          /* Receive the foreign parts to compute BH accretion rates and do the
           * swallowing */
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_rho);
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_part_swallow);
          /* scheduler_activate_recv(s, cj->mpi.recv,
           * task_subtype_bpart_merger); */

          /* Send the local BHs to tag the particles to swallow and to do
           * feedback */
          scheduler_activate_send(s, ci->mpi.send, task_subtype_bpart_rho,
                                  cj_nodeID);
          scheduler_activate_send(s, ci->mpi.send, task_subtype_bpart_feedback,
                                  cj_nodeID);
          /* scheduler_activate_send(s, ci->mpi.send, task_subtype_bpart_merger,
           */
          /*                         cj_nodeID); */

          /* Drift before you send */
          cell_activate_drift_bpart(ci, s);
        }

        if (cj_active) {

          /* Receive the foreign BHs to tag particles to swallow and for
           * feedback */
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_bpart_rho);
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_bpart_feedback);

          /* Send the local part information */
          scheduler_activate_send(s, ci->mpi.send, task_subtype_rho, cj_nodeID);
          scheduler_activate_send(s, ci->mpi.send, task_subtype_part_swallow,
                                  cj_nodeID);

          /* if (ci->black_holes.count > 0)  */
          /* scheduler_activate_send(s, ci->mpi.send, task_subtype_bpart_merger,
           */
          /* 			  ci_nodeID); */

          /* Drift the cell which will be sent; note that not all sent
             particles will be drifted, only those that are needed. */
          cell_activate_drift_part(ci, s);
          /* if (ci->black_holes.count > 0) cell_activate_drift_bpart(ci, s); */
        }
      }
#endif
    }
  }

  /* Un-skip the swallow tasks involved with this cell. */
  for (struct link *l = c->black_holes.swallow; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const int ci_active = cell_is_active_black_holes(ci, e);
    const int cj_active = (cj != NULL) ? cell_is_active_black_holes(cj, e) : 0;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active || cj_active) &&
        (ci_nodeID == nodeID || cj_nodeID == nodeID)) {

      scheduler_activate(s, t);
    }
  }

  /* Un-skip the swallow tasks involved with this cell. */
  for (struct link *l = c->black_holes.do_gas_swallow; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const int ci_active = cell_is_active_black_holes(ci, e);
    const int cj_active = (cj != NULL) ? cell_is_active_black_holes(cj, e) : 0;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active || cj_active) &&
        (ci_nodeID == nodeID || cj_nodeID == nodeID)) {

      scheduler_activate(s, t);
    }
  }

  /* Un-skip the swallow tasks involved with this cell. */
  for (struct link *l = c->black_holes.do_bh_swallow; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const int ci_active = cell_is_active_black_holes(ci, e);
    const int cj_active = (cj != NULL) ? cell_is_active_black_holes(cj, e) : 0;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active || cj_active) &&
        (ci_nodeID == nodeID || cj_nodeID == nodeID)) {

      scheduler_activate(s, t);
    }
  }

  /* Un-skip the feedback tasks involved with this cell. */
  for (struct link *l = c->black_holes.feedback; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
    const int ci_active = cell_is_active_black_holes(ci, e);
    const int cj_active = (cj != NULL) ? cell_is_active_black_holes(cj, e) : 0;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active || cj_active) &&
        (ci_nodeID == nodeID || cj_nodeID == nodeID)) {

      scheduler_activate(s, t);

      if (t->type == task_type_pair || t->type == task_type_sub_pair) {
        /* Activate bh_out for each cell that is part of
         * a pair/sub_pair task as to not miss any dependencies */
        if (ci_nodeID == nodeID)
          scheduler_activate(s, ci->hydro.super->black_holes.black_holes_out);
        if (cj_nodeID == nodeID)
          scheduler_activate(s, cj->hydro.super->black_holes.black_holes_out);
      }
    }
  }

  /* Unskip all the other task types. */
  if (c->nodeID == nodeID && cell_is_active_black_holes(c, e)) {

    /* If the cell doesn't have any pair/sub_pair type tasks,
     * then we haven't unskipped all the implicit tasks yet. */
    if (c->black_holes.density_ghost != NULL)
      scheduler_activate(s, c->black_holes.density_ghost);
    if (c->black_holes.swallow_ghost_1 != NULL)
      scheduler_activate(s, c->black_holes.swallow_ghost_1);
    if (c->black_holes.swallow_ghost_2 != NULL)
      scheduler_activate(s, c->black_holes.swallow_ghost_2);
    if (c->black_holes.swallow_ghost_3 != NULL)
      scheduler_activate(s, c->black_holes.swallow_ghost_3);
    // if (c->black_holes.black_holes_in != NULL)
    //   scheduler_activate(s, c->black_holes.black_holes_in);
    /* if (c->black_holes.black_holes_out != NULL) */
    /*   scheduler_activate(s, c->black_holes.black_holes_out); */
    if (c->kick1 != NULL) scheduler_activate(s, c->kick1);
    if (c->kick2 != NULL) scheduler_activate(s, c->kick2);
    if (c->timestep != NULL) scheduler_activate(s, c->timestep);
    if (c->top->timestep_collect != NULL)
      scheduler_activate(s, c->top->timestep_collect);
#ifdef WITH_CSDS
    if (c->csds != NULL) scheduler_activate(s, c->csds);
#endif
  }

  return rebuild;
}

/**
 * @brief Un-skips all the sinks tasks associated with a given cell and
 * checks if the space needs to be rebuilt.
 *
 * @param c the #cell.
 * @param s the #scheduler.
 *
 * @return 1 If the space needs rebuilding. 0 otherwise.
 */
int cell_unskip_sinks_tasks(struct cell *c, struct scheduler *s) {

  struct engine *e = s->space->e;
  const int with_timestep_sync = (e->policy & engine_policy_timestep_sync);
  const int nodeID = e->nodeID;
  int rebuild = 0;

  if (c->sinks.drift != NULL)
    if (cell_is_active_sinks(c, e) || cell_is_active_hydro(c, e)) {
      cell_activate_drift_sink(c, s);
    }

  /* Un-skip the star formation tasks involved with this cell. */
  for (struct link *l = c->sinks.swallow; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    const int ci_active =
        cell_is_active_sinks(ci, e) || cell_is_active_hydro(ci, e);

    const int cj_active = (cj != NULL) && (cell_is_active_sinks(cj, e) ||
                                           cell_is_active_hydro(cj, e));

    /* Activate the drifts */
    if (t->type == task_type_self && ci_active) {
      cell_activate_drift_sink(ci, s);
      cell_activate_drift_part(ci, s);
      cell_activate_sink_formation_tasks(ci->top, s);
      if (with_timestep_sync) cell_activate_sync_part(ci, s);
    }

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active || cj_active) &&
        (ci_nodeID == nodeID || cj_nodeID == nodeID)) {
      scheduler_activate(s, t);

      if (t->type == task_type_pair) {
        /* For the mergers */
        if (cj_nodeID == nodeID) {
          cell_activate_drift_sink(cj, s);
          cell_activate_sink_formation_tasks(cj->top, s);
        }
        if (ci_nodeID == nodeID) {
          cell_activate_drift_sink(ci, s);
          if (ci->top != cj->top) {
            cell_activate_sink_formation_tasks(ci->top, s);
          }
        }

        /* Do ci */
        if (ci_active) {
          /* hydro for cj */
          atomic_or(&cj->hydro.requires_sorts, 1 << t->flags);
          cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;

          /* Activate the drift tasks. */
          if (cj_nodeID == nodeID) cell_activate_drift_part(cj, s);
          if (cj_nodeID == nodeID && with_timestep_sync)
            cell_activate_sync_part(cj, s);

          /* Check the sorts and activate them if needed. */
          cell_activate_hydro_sorts(cj, t->flags, s);
        }

        /* Do cj */
        if (cj_active) {
          /* hydro for ci */
          atomic_or(&ci->hydro.requires_sorts, 1 << t->flags);
          ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;

          /* Activate the drift tasks. */
          if (ci_nodeID == nodeID) cell_activate_drift_part(ci, s);
          if (ci_nodeID == nodeID && with_timestep_sync)
            cell_activate_sync_part(ci, s);

          /* Check the sorts and activate them if needed. */
          cell_activate_hydro_sorts(ci, t->flags, s);
        }
      }

      else if (t->type == task_type_sub_self) {
        cell_activate_subcell_sinks_tasks(ci, NULL, s, with_timestep_sync);
      }

      else if (t->type == task_type_sub_pair) {
        cell_activate_subcell_sinks_tasks(ci, cj, s, with_timestep_sync);
      }
    }

    /* Only interested in pair interactions as of here. */
    if (t->type == task_type_pair || t->type == task_type_sub_pair) {
      /* Check whether there was too much particle motion, i.e. the
         cell neighbour conditions were violated. */
      if (cell_need_rebuild_for_sinks_pair(ci, cj)) rebuild = 1;
      if (cell_need_rebuild_for_sinks_pair(cj, ci)) rebuild = 1;

      /* Activate all sink_in tasks for each cell involved
       * in pair/sub_pair type tasks */
      if (ci_nodeID == nodeID)
        scheduler_activate(s, ci->hydro.super->sinks.sink_in);
      if (cj_nodeID == nodeID)
        scheduler_activate(s, cj->hydro.super->sinks.sink_in);

#ifdef WITH_MPI
      error("TODO");
#endif
    }
  }

  for (struct link *l = c->sinks.do_sink_swallow; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    const int ci_active =
        cell_is_active_sinks(ci, e) || cell_is_active_hydro(ci, e);

    const int cj_active = (cj != NULL) && (cell_is_active_sinks(cj, e) ||
                                           cell_is_active_hydro(cj, e));

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active || cj_active) &&
        (ci_nodeID == nodeID || cj_nodeID == nodeID)) {
      scheduler_activate(s, t);
    }
  }

  for (struct link *l = c->sinks.do_gas_swallow; l != NULL; l = l->next) {
    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    const int ci_active =
        cell_is_active_sinks(ci, e) || cell_is_active_hydro(ci, e);

    const int cj_active = (cj != NULL) && (cell_is_active_sinks(cj, e) ||
                                           cell_is_active_hydro(cj, e));

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active || cj_active) &&
        (ci_nodeID == nodeID || cj_nodeID == nodeID)) {
      scheduler_activate(s, t);

      if (t->type == task_type_pair || t->type == task_type_sub_pair) {
        /* Activate sinks_out for each cell that is part of
         * a pair/sub_pair task as to not miss any dependencies */
        if (ci_nodeID == nodeID)
          scheduler_activate(s, ci->hydro.super->sinks.sink_out);
        if (cj_nodeID == nodeID)
          scheduler_activate(s, cj->hydro.super->sinks.sink_out);
      }
    }
  }

  /* Unskip all the other task types. */
  if (c->nodeID == nodeID &&
      (cell_is_active_sinks(c, e) || cell_is_active_hydro(c, e))) {

    if (c->sinks.sink_in != NULL) scheduler_activate(s, c->sinks.sink_in);
    if (c->sinks.sink_ghost1 != NULL)
      scheduler_activate(s, c->sinks.sink_ghost1);
    if (c->sinks.sink_ghost2 != NULL)
      scheduler_activate(s, c->sinks.sink_ghost2);
    if (c->sinks.sink_out != NULL) scheduler_activate(s, c->sinks.sink_out);
    if (c->kick1 != NULL) scheduler_activate(s, c->kick1);
    if (c->kick2 != NULL) scheduler_activate(s, c->kick2);
    if (c->timestep != NULL) scheduler_activate(s, c->timestep);
    if (c->top->timestep_collect != NULL)
      scheduler_activate(s, c->top->timestep_collect);
#ifdef WITH_CSDS
    if (c->csds != NULL) scheduler_activate(s, c->csds);
#endif
  }

  return rebuild;
}

/**
 * @brief Un-skips all the RT tasks associated with a given cell and checks
 * if the space needs to be rebuilt.
 *
 * @param c the #cell.
 * @param s the #scheduler.
 * @param sub_cycle 1 if this is unskipping during an RT subcycle, 0 if normal
 * unskip
 *
 * @return 1 If the space needs rebuilding. 0 otherwise.
 */
int cell_unskip_rt_tasks(struct cell *c, struct scheduler *s,
                         const int sub_cycle) {

  /* Do we have work here? */
  if (c->hydro.count == 0) return 0;

  struct engine *e = s->space->e;
  const int nodeID = e->nodeID;
  int rebuild = 0; /* TODO: implement rebuild conditions? */

  /* Note: we only get this far if engine_policy_rt is flagged. */
  if (!(e->policy & engine_policy_rt)) error("Unskipping RT tasks without RT");

  for (struct link *l = c->rt.rt_gradient; l != NULL; l = l->next) {

    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif
    const int ci_active = cell_is_rt_active(ci, e);
    const int cj_active = (cj != NULL) && cell_is_rt_active(cj, e);

    /* Only activate tasks that involve a local active cell. */
    if ((ci_active && ci_nodeID == nodeID) ||
        (cj_active && cj_nodeID == nodeID)) {
      scheduler_activate(s, t);

      if (!sub_cycle) {
        /* Activate sorts only during main/normal steps. */
        if (t->type == task_type_pair) {
          atomic_or(&ci->hydro.requires_sorts, 1 << t->flags);
          atomic_or(&cj->hydro.requires_sorts, 1 << t->flags);
          ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;
          cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;

          /* Check the sorts and activate them if needed. */
          cell_activate_rt_sorts(ci, t->flags, s);
          cell_activate_rt_sorts(cj, t->flags, s);
        }

        /* Store current values of dx_max and h_max. */
        else if (t->type == task_type_sub_self) {
          cell_activate_subcell_rt_tasks(ci, NULL, s, sub_cycle);
        }

        /* Store current values of dx_max and h_max. */
        else if (t->type == task_type_sub_pair) {
          cell_activate_subcell_rt_tasks(ci, cj, s, sub_cycle);
        }
      }
    }

    /* Only interested in pair interactions as of here. */
    if (t->type == task_type_pair || t->type == task_type_sub_pair) {

#ifdef WITH_MPI

      /* Activate the send/recv tasks. */
      if (ci_nodeID != nodeID) {
        /* If the local cell is active, receive data from the foreign cell. */
        if (cj_active) {
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_rt_gradient);
          if (sub_cycle) {
            /* If we're in a sub-cycle, then there should be no sorts. But since
             * hydro sorts won't be active then, the RT sorts would run. Make
             * sure the cells are also marked to skip the RT sorts, otherwise
             * the 'sorted' flags will be wrongly set after a recv rt_gradient.
             * The recv tasks might also run on a higher level than the current
             * cell, so walk all the way up. */
            cell_set_skip_rt_sort_flag_up(ci);
          }

          /* We only need updates later on if the other cell is active too */
          if (ci_active) {
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_rt_transport);
          }
        } else if (ci_active) {
#ifdef MPI_SYMMETRIC_FORCE_INTERACTION
          /* If the local cell is inactive and the remote cell is active, we
           * still need to receive stuff to be able to do the force interaction
           * on this node as well.
           * The gradient recv is only necessary in normal steps in case we need
           * to sort, not during sub-cycles. */
          if (!sub_cycle)
            scheduler_activate_recv(s, ci->mpi.recv, task_subtype_rt_gradient);
          scheduler_activate_recv(s, ci->mpi.recv, task_subtype_rt_transport);
          if (sub_cycle) cell_set_skip_rt_sort_flag_up(ci);
#endif
        }

        /* Is the foreign cell active and will need stuff from us? */
        if (ci_active) {

          scheduler_activate_send(s, cj->mpi.send, task_subtype_rt_gradient,
                                  ci_nodeID);

          if (cj_active) {
            scheduler_activate_send(s, cj->mpi.send, task_subtype_rt_transport,
                                    ci_nodeID);
          }
        } else if (cj_active) {
#ifdef MPI_SYMMETRIC_FORCE_INTERACTION
          /* If the foreign cell is inactive, but the local cell is active,
           * we still need to send stuff to be able to do the force interaction
           * on both nodes.
           * The gradient send is only necessary in normal steps in case we need
           * to sort, not during sub-cycles. */
          if (!sub_cycle)
            scheduler_activate_send(s, cj->mpi.send, task_subtype_rt_gradient,
                                    ci_nodeID);
          scheduler_activate_send(s, cj->mpi.send, task_subtype_rt_transport,
                                  ci_nodeID);
#endif
        }

      } else if (cj_nodeID != nodeID) {

        /* If the local cell is active, receive data from the foreign cell. */
        if (ci_active) {
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_rt_gradient);
          if (sub_cycle) {
            /* No RT sorts during sub-cycling */
            cell_set_skip_rt_sort_flag_up(cj);
          }

          /* We only need updates later on if the other cell is active too */
          if (cj_active) {
            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_rt_transport);
          }
        } else if (cj_active) {
#ifdef MPI_SYMMETRIC_FORCE_INTERACTION
          /* If the local cell is inactive and the remote cell is active, we
           * still need to receive stuff to be able to do the force interaction
           * on this node as well.
           * The gradient recv is only necessary in normal steps in case we need
           * to sort, not during sub-cycles. */
          if (!sub_cycle)
            scheduler_activate_recv(s, cj->mpi.recv, task_subtype_rt_gradient);
          scheduler_activate_recv(s, cj->mpi.recv, task_subtype_rt_transport);
          if (sub_cycle) cell_set_skip_rt_sort_flag_up(cj);
#endif
        }

        /* Is the foreign cell active and will need stuff from us? */
        if (cj_active) {

          scheduler_activate_send(s, ci->mpi.send, task_subtype_rt_gradient,
                                  cj_nodeID);

          if (ci_active) {
            scheduler_activate_send(s, ci->mpi.send, task_subtype_rt_transport,
                                    cj_nodeID);
          }
        } else if (ci_active) {
#ifdef MPI_SYMMETRIC_FORCE_INTERACTION
          /* If the foreign cell is inactive, but the local cell is active,
           * we still need to send stuff to be able to do the force interaction
           * on both nodes
           * The gradient send is only necessary in normal steps in case we need
           * to sort, not during sub-cycles. */
          if (!sub_cycle)
            scheduler_activate_send(s, ci->mpi.send, task_subtype_rt_gradient,
                                    cj_nodeID);
          scheduler_activate_send(s, ci->mpi.send, task_subtype_rt_transport,
                                  cj_nodeID);
#endif
        }
      }
#endif
    }
  }

  for (struct link *l = c->rt.rt_transport; l != NULL; l = l->next) {

    struct task *t = l->t;
    struct cell *ci = t->ci;
    struct cell *cj = t->cj;
#ifdef WITH_MPI
    const int ci_nodeID = ci->nodeID;
    const int cj_nodeID = (cj != NULL) ? cj->nodeID : -1;
#else
    const int ci_nodeID = nodeID;
    const int cj_nodeID = nodeID;
#endif

    const int ci_active = cell_is_rt_active(ci, e);
    const int cj_active = ((cj != NULL) && cell_is_rt_active(cj, e));

    if ((ci_active && ci_nodeID == nodeID) ||
        (cj_active && cj_nodeID == nodeID)) {
      scheduler_activate(s, t);

      if (t->type == task_type_pair || t->type == task_type_sub_pair) {

        /* Activate transport_out for each cell that is part of
         * a pair/sub_pair task as to not miss any dependencies */
        if (ci_nodeID == nodeID)
          scheduler_activate(s, ci->hydro.super->rt.rt_transport_out);
        if (cj_nodeID == nodeID)
          scheduler_activate(s, cj->hydro.super->rt.rt_transport_out);
      }
    }
  }

  /* Unskip all the other task types */
  if (cell_is_rt_active(c, e)) {
    if (c->nodeID == nodeID) {
      if (c->rt.rt_in != NULL) scheduler_activate(s, c->rt.rt_in);
      if (c->rt.rt_ghost1 != NULL) scheduler_activate(s, c->rt.rt_ghost1);
      if (c->rt.rt_ghost2 != NULL) scheduler_activate(s, c->rt.rt_ghost2);
      if (c->rt.rt_transport_out != NULL)
        scheduler_activate(s, c->rt.rt_transport_out);
      if (c->rt.rt_tchem != NULL) scheduler_activate(s, c->rt.rt_tchem);
      if (c->rt.rt_out != NULL) scheduler_activate(s, c->rt.rt_out);
    } else {
#if defined(MPI_SYMMETRIC_FORCE_INTERACTION) && defined(WITH_MPI)
      /* Additionally unskip force interactions between inactive local cell and
       * active remote cell. (The cell unskip will only be called for active
       * cells, so, we have to do this now, from the active remote cell). */
      for (struct link *l = c->rt.rt_transport; l != NULL; l = l->next) {
        struct task *t = l->t;
        if (t->type != task_type_pair && t->type != task_type_sub_pair)
          continue;

        struct cell *ci = l->t->ci;
        struct cell *cj = l->t->cj;

        const int ci_active = cell_is_rt_active(ci, e);
        const int cj_active = cell_is_rt_active(cj, e);
        const int ci_nodeID = ci->nodeID;
        const int cj_nodeID = cj->nodeID;
        if ((!ci_active && ci_nodeID == nodeID && cj_active &&
             cj_nodeID != nodeID) ||
            (!cj_active && cj_nodeID == nodeID && ci_active &&
             ci_nodeID != nodeID)) {
          scheduler_activate(s, l->t);

          if (!sub_cycle) {
            /* Activate sorts only during main/normal steps. */
            if (t->type == task_type_pair) {
              atomic_or(&ci->hydro.requires_sorts, 1 << t->flags);
              atomic_or(&cj->hydro.requires_sorts, 1 << t->flags);
              ci->hydro.dx_max_sort_old = ci->hydro.dx_max_sort;
              cj->hydro.dx_max_sort_old = cj->hydro.dx_max_sort;

              /* Check the sorts and activate them if needed. */
              cell_activate_rt_sorts(ci, t->flags, s);
              cell_activate_rt_sorts(cj, t->flags, s);
            }

            /* Store current values of dx_max and h_max. */
            else if (t->type == task_type_sub_pair) {
              cell_activate_subcell_rt_tasks(ci, cj, s, sub_cycle);
            }
          }
        }
      }
#endif
    }

    /* The rt_advance_cell_time tasks also run on foreign cells */
    if (c->super != NULL && c->super->rt.rt_advance_cell_time != NULL) {
      scheduler_activate(s, c->super->rt.rt_advance_cell_time);
    }
    if (sub_cycle) {
      /* The rt_collect_times tasks replace the timestep_collect tasks
       * during sub-cycles, so we only activate it when sub-cycling. */
      if (c->top->rt.rt_collect_times != NULL)
        scheduler_activate(s, c->top->rt.rt_collect_times);
    } else {
      /* Otherwise, make sure timestep_collect and timestep is active. */
      if (c->top->timestep_collect != NULL)
        scheduler_activate(s, c->top->timestep_collect);
    }
  }

  return rebuild;
}
