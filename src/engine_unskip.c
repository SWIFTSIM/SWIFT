/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include "engine.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "memswap.h"

/* Load the profiler header, if needed. */
#ifdef WITH_PROFILER
#include <gperftools/profiler.h>
#endif

/**
 * @brief Broad categories of tasks.
 *
 * Each category is unskipped independently
 * of the others.
 */
enum task_broad_types {
  task_broad_types_hydro = 1,
  task_broad_types_gravity,
  task_broad_types_stars,
  task_broad_types_sinks,
  task_broad_types_black_holes,
  task_broad_types_rt,
  task_broad_types_count,
};

/**
 * @brief Meta-data for the unskipping
 */
struct unskip_data {

  /*! The #engine */
  struct engine *e;

  /*! Pointer to the start of the list of cells to unskip */
  int *list_base;

  /*! Number of times the list has been duplicated */
  int multiplier;

  /*! The number of active cells (without dulication) */
  int num_active_cells;

  /*! The #task_broad_types corresponding to each copy of the list */
  enum task_broad_types task_types[task_broad_types_count];
};

/**
 * @brief Unskip any hydro tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void engine_do_unskip_hydro(struct cell *c, struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure we only get zoom cells here when we are in zoom land. */
  if (e->s->with_zoom_region && c->type != cell_type_zoom) {
    error("Unskipping hydro tasks in a non-zoom cell (%s/%s) in zoom land.",
          cellID_names[c->type], subcellID_names[c->subtype]);
  }
#endif

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->hydro.count == 0) return;

#ifndef MPI_SYMMETRIC_FORCE_INTERACTION
  /* Skip inactive cells. */
  if (!cell_is_active_hydro(c, e)) return;
#endif

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_hydro(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_hydro_tasks(c, &e->sched);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any stars tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 * @param with_star_formation Are we running with star formation switched on?
 */
static void engine_do_unskip_stars(struct cell *c, struct engine *e,
                                   const int with_star_formation,
                                   const int with_star_formation_sink) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  const int non_empty =
      c->stars.count > 0 || (with_star_formation && c->hydro.count > 0) ||
      (with_star_formation_sink && (c->hydro.count > 0 || c->sinks.count > 0));

  /* Ignore empty cells. */
  if (!non_empty) return;

  const int ci_active = cell_need_activating_stars(c, e, with_star_formation,
                                                   with_star_formation_sink);

  /* Skip inactive cells. */
  if (!ci_active) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_stars(cp, e, with_star_formation,
                               with_star_formation_sink);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_stars_tasks(
      c, &e->sched, with_star_formation, with_star_formation_sink);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any black hole tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void engine_do_unskip_black_holes(struct cell *c, struct engine *e) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_black_holes(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_black_holes_tasks(c, &e->sched);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any sink tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void engine_do_unskip_sinks(struct cell *c, struct engine *e) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells. */
  if (c->sinks.count == 0 && c->hydro.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_sinks(c, e) && !cell_is_active_hydro(c, e)) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_sinks(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_sinks_tasks(c, &e->sched);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Is a non-empty top-level cell active for any of the broad task types
 * of the run?
 *
 * @param c The top-level #cell.
 * @param e The #engine.
 */
static int engine_unskip_top_cell_is_active(const struct cell *c,
                                            const struct engine *e) {

  const int nodeID = e->nodeID;
  const int with_hydro = e->policy & engine_policy_hydro;
  const int with_self_grav = e->policy & engine_policy_self_gravity;
  const int with_ext_grav = e->policy & engine_policy_external_gravity;
  const int with_stars = e->policy & engine_policy_stars;
  const int with_sinks = e->policy & engine_policy_sinks;
  const int with_feedback = e->policy & engine_policy_feedback;
  const int with_black_holes = e->policy & engine_policy_black_holes;
  const int with_rt = e->policy & engine_policy_rt;

  if (cell_is_empty(c)) return 0;

  return (with_hydro && cell_is_active_hydro(c, e)) ||
         (with_self_grav && cell_is_active_gravity(c, e)) ||
         (with_ext_grav && c->nodeID == nodeID &&
          cell_is_active_gravity(c, e)) ||
         (with_feedback && cell_is_active_stars(c, e)) ||
         (with_stars && c->nodeID == nodeID && cell_is_active_stars(c, e)) ||
         (with_sinks && cell_is_active_sinks(c, e)) ||
         (with_black_holes && cell_is_active_black_holes(c, e)) ||
         (with_rt && cell_is_rt_active(c, e));
}

/**
 * @brief Does the gravity unskip of a cell act on it?
 *
 * @param c The #cell.
 * @param e The #engine.
 */
static int engine_unskip_gravity_acts_on(const struct cell *c,
                                         const struct engine *e) {

  /* Below the level where tasks are, empty (but not void) or inactive. */
  return cell_get_flag(c, cell_flag_has_tasks) && !cell_is_empty_mpole(c) &&
         cell_is_active_gravity(c, e);
}

/**
 * @brief Does the gravity walk of a cell it acts on go down to its progeny?
 *
 * @param c The #cell.
 */
static int engine_unskip_gravity_enters_progeny(const struct cell *c) {

  return cell_is_split_or_void(c) && (cell_is_above_diff_grav_depth(c) ||
                                      c->grav.tasks_below_diff_grav_depth ||
                                      c->subtype == cell_subtype_void);
}

/**
 * @brief Does the gravity unskip of a cell recurse into its progeny?
 *
 * @param c The #cell.
 * @param e The #engine.
 */
static int engine_unskip_gravity_recurses(struct cell *c,
                                          const struct engine *e) {

  return engine_unskip_gravity_acts_on(c, e) &&
         engine_unskip_gravity_enters_progeny(c);
}

/**
 * @brief Unskip any gravity tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 * @param walk_below Also unskip the progeny of @p c?
 */
static void engine_do_unskip_gravity(struct cell *c, struct engine *e,
                                     const int walk_below) {

  if (!engine_unskip_gravity_acts_on(c, e)) return;

  /* At the top level we need to recursively check particles haven't moved
   * too far for the mesh gravity (if using the mesh). */
  if (c->depth == 0 && e->s->periodic) {

    /* In zoom land we need to redirect to the zoom version. */
    if (!e->s->with_zoom_region) {
      cell_check_grav_mesh_pairs(c, e);
    } else {
      cell_check_grav_mesh_pairs_zoom(c, e);
    }
  }

  /* Recurse */
  if (walk_below && engine_unskip_gravity_enters_progeny(c)) {
    for (int k = 0; k < 8; k++) {
      struct cell *cp = c->progeny[k];
      if (cp != NULL) engine_do_unskip_gravity(cp, e, /*walk_below=*/1);
    }
  }

  /* Unskip any active tasks. */
  cell_unskip_gravity_tasks(c, &e->sched);
}

/**
 * @brief Unskip any radiative transfer tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 * @param sub_cycle 1 if this is a call for a sub cycle, 0 otherwise
 */
static void engine_do_unskip_rt(struct cell *c, struct engine *e,
                                const int sub_cycle) {

  /* Note: we only get this far if engine_policy_rt is flagged. */
#ifdef SWIFT_DEBUG_CHECKS
  if (!(e->policy & engine_policy_rt))
    error("Unksipping RT stuff without the policy being on");
#endif

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Do we have work to do? */
  if (c->hydro.count == 0) return;
  if (!cell_is_rt_active(c, e)) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        engine_do_unskip_rt(c->progeny[k], e, sub_cycle);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_rt_tasks(c, &e->sched, sub_cycle);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/*! Number of tree levels by which a local top-level cell is split into jobs.
 * 2 gives at most 1 + 8 + 64 jobs per top-level cell. */
static const int engine_unskip_split_levels = 2;

/*! Chunk size of the gravity items. Chunk 1 costs one atomic per item and
 * automatic chunking gives one thread a long run of the heavy items of a
 * single top-level cell, which are neighbours in the list. */
static const int engine_unskip_gravity_chunk = 8;

/**
 * @brief A cell to unskip from one job.
 */
struct unskip_item {

  /*! The #cell */
  struct cell *c;

  /*! Does the job also unskip the progeny of the cell? */
  int walk_below;
};

/**
 * @brief A growable list of #unskip_item and the callbacks that build it.
 */
struct unskip_item_list {

  /*! The #engine */
  struct engine *e;

  /*! Does the unskip act on a cell at all? */
  int (*acts_on)(const struct cell *c, const struct engine *e);

  /*! Does the walk of a cell recurse into its progeny? */
  int (*recurses)(struct cell *c, const struct engine *e);

  /*! Is a progeny unskipped by a job of its own? Can be NULL. */
  int (*has_own_job)(const struct cell *c, struct cell *cp,
                     const struct engine *e);

  /*! The per-cell unskip action of a job */
  void (*unskip)(struct cell *c, struct engine *e, int walk_below);

  /*! The items */
  struct unskip_item *items;

  /*! Number of items and allocated size */
  int count, size;
};

/**
 * @brief Append an item to an #unskip_item_list.
 *
 * @param list The #unskip_item_list.
 * @param c The #cell of the item.
 * @param walk_below Does the item also unskip the progeny of @p c?
 */
static void engine_unskip_add_item(struct unskip_item_list *list,
                                   struct cell *c, const int walk_below) {

  if (list->count == list->size) {
    list->size = list->size > 0 ? 2 * list->size : 1024;
    list->items = (struct unskip_item *)realloc(
        list->items, list->size * sizeof(struct unskip_item));
    if (list->items == NULL) error("Couldn't allocate the unskip items.");
  }
  list->items[list->count].c = c;
  list->items[list->count].walk_below = walk_below;
  list->count++;
}

/**
 * @brief Split the walk of a cell into items, down to a fixed number of levels.
 *
 * A local cell whose walk recurses becomes an item that unskips the cell alone
 * and its progeny are split in turn. The other cells become items that walk
 * their whole subtree.
 *
 * @param list The #unskip_item_list to fill.
 * @param c The #cell.
 * @param levels The number of levels left to split.
 */
static void engine_unskip_split_cell(struct unskip_item_list *list,
                                     struct cell *c, const int levels) {

  const struct engine *e = list->e;

  if (!list->acts_on(c, e)) return;

  if (levels == 0 || c->nodeID != e->nodeID || !list->recurses(c, e)) {
    engine_unskip_add_item(list, c, /*walk_below=*/1);
    return;
  }

  engine_unskip_add_item(list, c, /*walk_below=*/0);
  for (int k = 0; k < 8; k++) {
    struct cell *cp = c->progeny[k];
    if (cp == NULL) continue;
    if (list->has_own_job != NULL && list->has_own_job(c, cp, e)) continue;
    engine_unskip_split_cell(list, cp, levels - 1);
  }
}

/**
 * @brief Build the items of a set of top-level cells.
 *
 * @param list The #unskip_item_list with its callbacks set and no items.
 * @param cell_ids The indices of the top-level cells.
 * @param num_cells The number of top-level cells.
 */
static void engine_unskip_build_items(struct unskip_item_list *list,
                                      const int *cell_ids,
                                      const int num_cells) {

  struct cell *cells_top = list->e->s->cells_top;
  for (int k = 0; k < num_cells; k++)
    engine_unskip_split_cell(list, &cells_top[cell_ids[k]],
                             engine_unskip_split_levels);
}

/**
 * @brief Mapper function to run the unskip action of a list of items.
 *
 * @param map_data An array of #unskip_item.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to the #unskip_item_list.
 */
static void engine_do_unskip_items_mapper(void *map_data, int num_elements,
                                          void *extra_data) {

  const struct unskip_item_list *list =
      (const struct unskip_item_list *)extra_data;
  const struct unskip_item *items = (const struct unskip_item *)map_data;

  for (int ind = 0; ind < num_elements; ind++)
    list->unskip(items[ind].c, list->e, items[ind].walk_below);
}

/**
 * @brief Mapper function to unskip active tasks.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an unskip_data structure.
 */
void engine_do_unskip_mapper(void *map_data, int num_elements,
                             void *extra_data) {

  /* Unpack the meta data */
  struct unskip_data *data = (struct unskip_data *)extra_data;
  const int num_active_cells = data->num_active_cells;
  const enum task_broad_types *const task_types = data->task_types;
  const int *const list_base = data->list_base;
  struct engine *e = data->e;
  struct cell *const cells_top = e->s->cells_top;

  /* What policies are we running? */
  const int with_sinks = e->policy & engine_policy_sinks;
  const int with_stars = e->policy & engine_policy_stars;
  const int with_star_formation = e->policy & engine_policy_star_formation;
  const int with_star_formation_sink = with_sinks && with_stars;

  /* The current chunk of active cells */
  const int *const local_cells = (int *)map_data;

  /* Loop over this thread's chunk of cells to unskip */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Handle on the cell */
    struct cell *const c = &cells_top[local_cells[ind]];

    /* In what copy of the global list are we?
     * This gives us the broad type of task we are working on. */
    const ptrdiff_t delta = &local_cells[ind] - list_base;
    const int type = delta / num_active_cells;

    /* When running a zoom background cells only have gravity tasks. Skip
     * here to enforce this. */
    if (c->type == cell_type_bkg &&
        task_types[type] != task_broad_types_gravity)
      continue;

#ifdef SWIFT_DEBUG_CHECKS
    if (type >= data->multiplier) error("Invalid broad task type!");
    if (c == NULL) error("Got an invalid cell index!");
#endif

    /* What broad type of tasks are we unskipping? */
    switch (task_types[type]) {
      case task_broad_types_hydro:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_hydro))
          error("Trying to unskip hydro tasks in a non-hydro run!");
#endif
        engine_do_unskip_hydro(c, e);
        break;
      case task_broad_types_stars:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_stars))
          error("Trying to unskip star tasks in a non-stars run!");
#endif
        engine_do_unskip_stars(c, e, with_star_formation,
                               with_star_formation_sink);
        break;
      case task_broad_types_sinks:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_sinks))
          error("Trying to unskip sink tasks in a non-sinks run!");
#endif
        engine_do_unskip_sinks(c, e);
        break;
      case task_broad_types_black_holes:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_black_holes))
          error("Trying to unskip black holes tasks in a non-BH run!");
#endif
        engine_do_unskip_black_holes(c, e);
        break;
      case task_broad_types_rt:
#ifdef SWIFT_DEBUG_CHECKS
        if (!(e->policy & engine_policy_rt))
          error("Trying to unskip radiative transfer tasks in a non-rt run!");
#endif
        engine_do_unskip_rt(c, e, /*sub_cycle=*/0);
        break;
      default:
#ifdef SWIFT_DEBUG_CHECKS
        error("Invalid broad task type!");
#endif
        continue;
    }
  }
}

#ifdef SWIFT_DEBUG_CHECKS
/**
 * @brief Does a top-level cell get its own job in engine_unskip()?
 *
 * This is the membership rule of space_list_useful_top_level_cells() combined
 * with the activity filter of engine_unskip().
 *
 * @param c The top-level #cell.
 * @param e The #engine.
 */
static int engine_unskip_top_cell_is_listed(struct cell *c,
                                            const struct engine *e) {

  return space_top_cell_has_tasks(c) && engine_unskip_top_cell_is_active(c, e);
}

/**
 * @brief Check that engine_unskip_top_cell_is_listed() holds exactly for the
 * top-level cells of the active list.
 *
 * @param e The #engine.
 * @param active_cells The indices of the active top-level cells.
 * @param num_active_cells The number of active top-level cells.
 */
static void engine_unskip_check_listed_rule(struct engine *e,
                                            const int *active_cells,
                                            const int num_active_cells) {

  struct space *s = e->s;
  char *listed = (char *)calloc(s->nr_cells, sizeof(char));
  if (listed == NULL) error("Couldn't allocate the active list flags.");

  for (int k = 0; k < num_active_cells; k++) listed[active_cells[k]] = 1;

  for (int cid = 0; cid < s->nr_cells; cid++) {
    struct cell *c = &s->cells_top[cid];
    if (engine_unskip_top_cell_is_listed(c, e) != listed[cid])
      error(
          "Active list rule does not match the active list (cellID=%lld "
          "listed=%d).",
          c->cellID, listed[cid]);
  }

  free(listed);
}

/*! Number of steps between two checks of the gravity unskip */
#ifndef ENGINE_UNSKIP_CHECK_PERIOD
#define ENGINE_UNSKIP_CHECK_PERIOD 1
#endif

/**
 * @brief Reference gravity unskip: one recursive walk per top-level cell.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void engine_do_unskip_gravity_reference(struct cell *c,
                                               struct engine *e) {

  /* Early abort (are we below the level where tasks are)? */
  if (!cell_get_flag(c, cell_flag_has_tasks)) return;

  /* Ignore empty cells but not void cells (in zoom land). */
  if (cell_is_empty_mpole(c)) return;

  /* Skip inactive cells. */
  if (!cell_is_active_gravity(c, e)) return;

  /* At the top level we need to recursively check particles haven't moved
   * too far for the mesh gravity (if using the mesh). */
  if (c->depth == 0 && e->s->periodic) {

    /* In zoom land we need to redirect to the zoom version. */
    if (!e->s->with_zoom_region) {
      cell_check_grav_mesh_pairs(c, e);
    } else {
      cell_check_grav_mesh_pairs_zoom(c, e);
    }
  }

  /* Recurse */
  if (cell_is_split_or_void(c) && (cell_is_above_diff_grav_depth(c) ||
                                   c->grav.tasks_below_diff_grav_depth ||
                                   c->subtype == cell_subtype_void)) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        engine_do_unskip_gravity_reference(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  cell_unskip_gravity_tasks(c, &e->sched);
}

/**
 * @brief What engine_unskip_check_flags() does on each cell.
 */
enum unskip_check_flags_mode {
  unskip_check_flags_count,
  unskip_check_flags_store,
  unskip_check_flags_compare,
  unskip_check_flags_clear,
  unskip_check_flags_restore,
};

/**
 * @brief Visit a cell and its progeny.
 *
 * @param c The #cell.
 * @param flags The stored flags, in the order of the visit.
 * @param ind The index of @p c in @p flags.
 * @param mode The #unskip_check_flags_mode.
 *
 * @return The index after the last cell visited.
 */
static long engine_unskip_check_flags_recursive(
    struct cell *c, uint32_t *flags, long ind,
    const enum unskip_check_flags_mode mode) {

  const uint32_t processed = cell_flag_unskip_self_grav_processed |
                             cell_flag_unskip_pair_grav_processed;
  const uint32_t drift = cell_flag_do_grav_drift | cell_flag_do_grav_sub_drift;
  const uint32_t reset = processed | drift;

  switch (mode) {
    case unskip_check_flags_count:
      break;
    case unskip_check_flags_store:
      flags[ind] = c->flags;
      break;
    case unskip_check_flags_compare: {
      const uint32_t added = c->flags & ~flags[ind] & ~reset;
      if (added)
        error(
            "The reference gravity unskip sets new cell flags 0x%x "
            "(cellID=%lld depth=%d nodeID=%d).",
            added, c->cellID, c->depth, c->nodeID);
      const uint32_t drift_diff = (c->flags ^ flags[ind]) & drift;
      if (drift_diff)
        error(
            "The reference gravity unskip and the split one set different "
            "drift flags 0x%x (cellID=%lld depth=%d nodeID=%d).",
            drift_diff, c->cellID, c->depth, c->nodeID);
      break;
    }
    case unskip_check_flags_clear:
      c->flags &= ~reset;
      break;
    case unskip_check_flags_restore:
      c->flags = (c->flags & ~reset) | (flags[ind] & reset);
      break;
  }
  ind++;

  for (int k = 0; k < 8; k++) {
    struct cell *cp = c->progeny[k];
    if (cp == NULL) continue;
    ind = engine_unskip_check_flags_recursive(cp, flags, ind, mode);
  }
  return ind;
}

/**
 * @brief Visit all the cells, see engine_unskip_check_flags_recursive().
 *
 * @param s The #space.
 * @param flags The stored flags, in the order of the visit.
 * @param mode The #unskip_check_flags_mode.
 *
 * @return The number of cells visited.
 */
static long engine_unskip_check_flags(struct space *s, uint32_t *flags,
                                      const enum unskip_check_flags_mode mode) {

  long ind = 0;
  for (int cid = 0; cid < s->nr_cells; cid++)
    ind = engine_unskip_check_flags_recursive(&s->cells_top[cid], flags, ind,
                                              mode);
  return ind;
}

/**
 * @brief Check that the reference gravity unskip activates no task and sets
 * no cell flag that the split gravity unskip did not.
 *
 * The processed and drift flags of the split unskip are cleared before the
 * reference walk, so that it walks the whole tree and sets the drift flags
 * itself, and restored after. The reference walk must then set the same drift
 * flags as the split one.
 *
 * @param e The #engine.
 * @param active_cells The indices of the active top-level cells.
 * @param num_active_cells The number of active top-level cells.
 */
static void engine_unskip_check_gravity(struct engine *e,
                                        const int *active_cells,
                                        const int num_active_cells) {

  struct space *s = e->s;
  struct scheduler *sched = &e->sched;

  const long num_cells =
      engine_unskip_check_flags(s, NULL, unskip_check_flags_count);
  uint32_t *flags = (uint32_t *)malloc(num_cells * sizeof(uint32_t));
  if (flags == NULL) error("Couldn't allocate the cell flags copy.");
  engine_unskip_check_flags(s, flags, unskip_check_flags_store);
  engine_unskip_check_flags(s, flags, unskip_check_flags_clear);

  const int active_count = sched->active_count;
  const int forcerebuild = e->forcerebuild;

  for (int k = 0; k < num_active_cells; k++)
    engine_do_unskip_gravity_reference(&s->cells_top[active_cells[k]], e);

  if (sched->active_count != active_count) {
    const struct task *t = &sched->tasks[sched->tid_active[active_count]];
    error(
        "The reference gravity unskip activates %d more tasks, first: task "
        "%lld %s/%s ci cellID=%lld cj cellID=%lld.",
        sched->active_count - active_count, (long long)(t - sched->tasks),
        taskID_names[t->type], subtaskID_names[t->subtype],
        t->ci != NULL ? t->ci->cellID : -1LL,
        t->cj != NULL ? t->cj->cellID : -1LL);
  }
  if (forcerebuild == 0 && e->forcerebuild != 0)
    error("The reference gravity unskip triggers a rebuild.");
  e->forcerebuild = forcerebuild;

  engine_unskip_check_flags(s, flags, unskip_check_flags_compare);
  engine_unskip_check_flags(s, flags, unskip_check_flags_restore);

  free(flags);
}
#endif

/**
 * @brief Unskip all the tasks that act on active cells at this time.
 *
 * @param e The #engine.
 */
void engine_unskip(struct engine *e) {

  const ticks tic = getticks();
  struct space *s = e->s;

  const int with_hydro = e->policy & engine_policy_hydro;
  const int with_self_grav = e->policy & engine_policy_self_gravity;
  const int with_ext_grav = e->policy & engine_policy_external_gravity;
  const int with_stars = e->policy & engine_policy_stars;
  const int with_sinks = e->policy & engine_policy_sinks;
  const int with_feedback = e->policy & engine_policy_feedback;
  const int with_black_holes = e->policy & engine_policy_black_holes;
  const int with_rt = e->policy & engine_policy_rt;

#ifdef WITH_PROFILER
  static int count = 0;
  char filename[100];
  sprintf(filename, "/tmp/swift_engine_do_usnkip_mapper_%06i.prof", count++);
  ProfilerStart(filename);
#endif  // WITH_PROFILER

  /* Move the active local cells to the top of the list. */
  int *local_cells = e->s->local_cells_with_tasks_top;
  int num_active_cells = 0;
  for (int k = 0; k < s->nr_local_cells_with_tasks; k++) {
    struct cell *c = &s->cells_top[local_cells[k]];

    if (engine_unskip_top_cell_is_active(c, e)) {

      if (num_active_cells != k)
        memswap(&local_cells[k], &local_cells[num_active_cells], sizeof(int));
      num_active_cells += 1;
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  engine_unskip_check_listed_rule(e, local_cells, num_active_cells);
#endif

  /* What kind of tasks do we have? */
  struct unskip_data data;
  bzero(&data, sizeof(struct unskip_data));
  int multiplier = 0;
  if (with_hydro) {
    data.task_types[multiplier] = task_broad_types_hydro;
    multiplier++;
  }
  if (with_feedback || with_stars) {
    data.task_types[multiplier] = task_broad_types_stars;
    multiplier++;
  }
  if (with_sinks) {
    data.task_types[multiplier] = task_broad_types_sinks;
    multiplier++;
  }
  if (with_black_holes) {
    data.task_types[multiplier] = task_broad_types_black_holes;
    multiplier++;
  }
  if (with_rt) {
    data.task_types[multiplier] = task_broad_types_rt;
    multiplier++;
  }

  /* Should we duplicate the list of active cells to better parallelise the
     unskip over the threads ? */
  int *local_active_cells;
  if (multiplier > 1) {

    /* Make space for copies of the list */
    local_active_cells =
        (int *)malloc(multiplier * num_active_cells * sizeof(int));
    if (local_active_cells == NULL)
      error(
          "Couldn't allocate memory for duplicated list of local active "
          "cells.");

    /* Make blind copies of the list */
    for (int m = 0; m < multiplier; m++) {
      memcpy(local_active_cells + m * num_active_cells, local_cells,
             num_active_cells * sizeof(int));
    }
  } else {
    local_active_cells = local_cells;
  }

  /* We now have a list of local active cells duplicated as many times as
   * we have broad task types. We can now release all the threads on the list */

  data.e = e;
  data.list_base = local_active_cells;
  data.num_active_cells = num_active_cells;
  data.multiplier = multiplier;

  /* Activate all the regular tasks (gravity-only runs have none) */
  if (multiplier > 0)
    threadpool_map(&e->threadpool, engine_do_unskip_mapper, local_active_cells,
                   num_active_cells * multiplier, sizeof(int), /*chunk=*/1,
                   &data);

  /* The gravity walks of large top-level cells are split into several jobs. */
  if (with_self_grav || with_ext_grav) {
    struct unskip_item_list grav_items;
    bzero(&grav_items, sizeof(struct unskip_item_list));
    grav_items.e = e;
    grav_items.acts_on = engine_unskip_gravity_acts_on;
    grav_items.recurses = engine_unskip_gravity_recurses;
    grav_items.has_own_job = NULL;
    grav_items.unskip = engine_do_unskip_gravity;
    const ticks tic_items = getticks();
    engine_unskip_build_items(&grav_items, local_cells, num_active_cells);
    const ticks toc_items = getticks();

    threadpool_map(&e->threadpool, engine_do_unskip_items_mapper,
                   grav_items.items, grav_items.count,
                   sizeof(struct unskip_item), engine_unskip_gravity_chunk,
                   &grav_items);

#ifdef SWIFT_DEBUG_CHECKS
    if (e->step % ENGINE_UNSKIP_CHECK_PERIOD == 0)
      engine_unskip_check_gravity(e, local_cells, num_active_cells);
#endif

    if (e->verbose)
      message(
          "gravity: %d active top-level cells split into %d jobs in %.3f %s.",
          num_active_cells, grav_items.count,
          clocks_from_ticks(toc_items - tic_items), clocks_getunit());

    free(grav_items.items);
  }

#ifdef WITH_PROFILER
  ProfilerStop();
#endif  // WITH_PROFILER

  /* Free stuff? */
  if (multiplier > 1) {
    free(local_active_cells);
  }

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void engine_do_unskip_sub_cycle_mapper(void *map_data, int num_elements,
                                       void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct cell *const cells_top = e->s->cells_top;

  /* The current chunk of active cells */
  const int *const local_cells = (int *)map_data;

  /* Loop over this thread's chunk of cells to unskip */
  for (int ind = 0; ind < num_elements; ind++) {

    /* Handle on the cell */
    struct cell *const c = &cells_top[local_cells[ind]];

    engine_do_unskip_rt(c, e, /*sub_cycle=*/1);
  }
}

/**
 * @brief Unskip all the RT tasks that are active during this sub-cycle.
 *
 * @param e The #engine.
 */
void engine_unskip_rt_sub_cycle(struct engine *e) {

  const ticks tic = getticks();
  struct space *s = e->s;
  const int with_rt = e->policy & engine_policy_rt;

  if (!with_rt) error("Unskipping sub-cycles when running without RT!");

  int *local_cells = e->s->local_cells_with_tasks_top;
  int num_active_cells = 0;
  for (int k = 0; k < s->nr_local_cells_with_tasks; k++) {
    struct cell *c = &s->cells_top[local_cells[k]];

    if (c->hydro.count == 0) continue;

    if (cell_is_rt_active(c, e)) {

      if (num_active_cells != k)
        memswap(&local_cells[k], &local_cells[num_active_cells], sizeof(int));
      num_active_cells += 1;
    }
  }

  /* Activate all the regular tasks */
  threadpool_map(&e->threadpool, engine_do_unskip_sub_cycle_mapper, local_cells,
                 num_active_cells, sizeof(int), /*chunk=*/1, e);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
