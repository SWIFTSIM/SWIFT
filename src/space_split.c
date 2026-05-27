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
#include "space.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "debug.h"
#include "engine.h"
#include "multipole.h"
#include "star_formation_logger.h"
#include "threadpool.h"

/* Depth budget per frontier entry. Each cell in the BFS frontier processes
 * this many levels depth-first before its surviving children re-enter the
 * frontier.  A value of 1 gives pure BFS; larger values trade dispatch
 * overhead for cache locality while preserving BFS load balance across
 * frontier cells. */
static const int dfs_levels_at_a_time = 3;

/* ------------------------------------------------------------------------- */
/* BFS-frontier infrastructure for space_split.                              */
/* ------------------------------------------------------------------------- */

/* Forward declarations of helpers defined further down this file. */
__attribute__((always_inline)) static INLINE int space_split_should_split(
    const struct space *s, const struct cell *c);
static void space_split_init_progeny(struct space *s, struct cell *c,
                                     const short int tpid);
static void space_split_finalize_leaf(struct space *s, struct cell *c);
static void space_split_fill_top_buffers(const struct cell *c,
                                         struct cell_buff *buff,
                                         struct cell_buff *sbuff,
                                         struct cell_buff *bbuff,
                                         struct cell_buff *gbuff,
                                         struct cell_buff *sink_buff);

/**
 * @brief A frontier of cells waiting to be split.
 *
 * Used as a flat array of @c cell pointers shared by all workers. Workers
 * read from @c cells and append produced children to a separate @e next
 * frontier whose @c count is grown atomically. Capacity is sized
 * conservatively to @c 8 * max_input_count before each level so no
 * reallocation is needed during the parallel append.
 */
struct space_split_frontier {

  /*! Pointers to the cells currently at this BFS level. */
  struct cell **cells;

  /*! Number of valid entries in @c cells. */
  int count;

  /*! Number of slots allocated in @c cells. */
  int capacity;
};

/**
 * @brief Reserve a slot in the next-level frontier and append a child
 *        cell pointer to it.
 *
 * Reservation is done with a single relaxed atomic increment of the
 * frontier @c count, which is safe because the buffer was sized to
 * @c 8 * (current frontier size) before the parallel append began, so
 * overflow is impossible.
 *
 * @param frontier The frontier being produced.
 * @param c The cell pointer to append.
 */
__attribute__((always_inline)) static INLINE void space_split_frontier_append(
    struct space_split_frontier *frontier, struct cell *c) {

  const int slot = atomic_inc(&frontier->count);
#ifdef SWIFT_DEBUG_CHECKS
  if (slot >= frontier->capacity)
    error("BFS frontier overflow: slot=%d capacity=%d", slot,
          frontier->capacity);
#endif
  frontier->cells[slot] = c;
}

/**
 * @brief Grow a frontier so that it can hold at least @c needed pointers.
 *
 * Used to pre-size the next-level frontier to @c 8 * size_of_current before a
 * parallel level so that #space_split_frontier_append needs no locking.
 * Frontiers are scratch storage that gets repopulated every level, so their
 * previous contents never need to be preserved when they grow.
 *
 * @param frontier The frontier to resize.
 * @param needed Required capacity (in pointers).
 */
static void space_split_frontier_ensure_capacity(
    struct space_split_frontier *frontier, const int needed) {

  if (needed <= frontier->capacity) return;
  struct cell **new_cells = NULL;
  if (swift_memalign("split_frontier", (void **)&new_cells,
                     SWIFT_STRUCT_ALIGNMENT,
                     sizeof(struct cell *) * needed) != 0)
    error("Failed to (re)allocate BFS frontier of size %d", needed);
  if (frontier->cells != NULL) swift_free("split_frontier", frontier->cells);
  frontier->cells = new_cells;
  frontier->capacity = needed;
}

/**
 * @brief Process a cell from the BFS frontier, recursing up to
 *        @c depth_remaining levels depth-first before returning survivors
 *        to the next frontier.
 *
 * Decides whether @c c should be split. If it stays a leaf,
 * #space_split_finalize_leaf is called to populate its summary fields. If
 * it splits, the 8 progeny are allocated, the cell's particles are
 * partitioned with @c cell_split, empty progeny are recycled, and each
 * surviving child is either recursed into (if @c depth_remaining > 0) or
 * atomically appended to @c next (if @c depth_remaining == 0).
 *
 * The caller passes the maximum number of BFS levels each frontier entry
 * should advance; at the deepest allowed level survivors become frontier
 * entries themselves. This gives a tunable balance between pure BFS
 * (@c depth_remaining == 0) and pure DFS (very large
 * @c depth_remaining), exploiting cache locality within a subtree while
 * still re-distributing work at the frontier level.
 *
 * @param s The #space.
 * @param top_buff, top_sbuff, top_bbuff, top_gbuff, top_sink_buff Transient
 *        per-top-cell split buffers (side arrays indexed by top-cell offset).
 * @param next The frontier being produced for the next BFS level.
 * @param c The cell to process.
 * @param tpid The threadpool tid of the calling worker.
 * @param depth_remaining Number of additional DFS levels to descend before
 *        enqueuing survivors (@c dfs_levels_at_a_time - 1 on entry from
 *        the mapper).
 */
static void space_split_process_subtree(
    struct space *s, struct cell_buff **top_buff,
    struct cell_buff **top_sbuff, struct cell_buff **top_bbuff,
    struct cell_buff **top_gbuff, struct cell_buff **top_sink_buff,
    struct space_split_frontier *next, struct cell *c, const short int tpid,
    const int depth_remaining) {

  /* Top-level cells inherit the worker's tpid (preserves the legacy
   * behaviour where a top cell's progeny share its tpid). */
  if (c->depth == 0) c->tpid = tpid;

  if (c->depth > space_cell_maxdepth) {
    error(
        "Exceeded maximum depth (%d) when splitting cells, aborting. This is "
        "most likely due to having too many particles at the exact same "
        "position, making the construction of a tree impossible.",
        space_cell_maxdepth);
  }

  /* Leaf? Finalise it and we're done with this cell. */
  if (!space_split_should_split(s, c)) {
    space_split_finalize_leaf(s, c);
    return;
  }

  /* Allocate and fill the transient split buffers lazily, only when a top
   * cell actually decides to split. */
  if (c->top == c) {
    const ptrdiff_t top_index = c - s->cells_top;

    if (top_buff[top_index] != NULL || top_sbuff[top_index] != NULL ||
        top_bbuff[top_index] != NULL || top_gbuff[top_index] != NULL ||
        top_sink_buff[top_index] != NULL)
      error("Top-level split buffers already allocated.");

    if (c->hydro.count > 0 &&
        swift_memalign("top_split_buff", (void **)&top_buff[top_index],
                       SWIFT_STRUCT_ALIGNMENT,
                       sizeof(struct cell_buff) * c->hydro.count) != 0)
      error("Failed to allocate hydro top-level split buffer.");

    if (c->stars.count > 0 &&
        swift_memalign("top_split_sbuff", (void **)&top_sbuff[top_index],
                       SWIFT_STRUCT_ALIGNMENT,
                       sizeof(struct cell_buff) * c->stars.count) != 0)
      error("Failed to allocate stars top-level split buffer.");

    if (c->black_holes.count > 0 &&
        swift_memalign("top_split_bbuff", (void **)&top_bbuff[top_index],
                       SWIFT_STRUCT_ALIGNMENT,
                       sizeof(struct cell_buff) * c->black_holes.count) != 0)
      error("Failed to allocate black-hole top-level split buffer.");

    if (c->grav.count > 0 &&
        swift_memalign("top_split_gbuff", (void **)&top_gbuff[top_index],
                       SWIFT_STRUCT_ALIGNMENT,
                       sizeof(struct cell_buff) * c->grav.count) != 0)
      error("Failed to allocate gravity top-level split buffer.");

    if (c->sinks.count > 0 &&
        swift_memalign("top_split_sink_buff",
                       (void **)&top_sink_buff[top_index],
                       SWIFT_STRUCT_ALIGNMENT,
                       sizeof(struct cell_buff) * c->sinks.count) != 0)
      error("Failed to allocate sink top-level split buffer.");

    space_split_fill_top_buffers(c, top_buff[top_index], top_sbuff[top_index],
                                 top_bbuff[top_index], top_gbuff[top_index],
                                 top_sink_buff[top_index]);
  }

  space_split_init_progeny(s, c, tpid);

  const ptrdiff_t top_index = c->top - s->cells_top;
  struct cell_buff *buff =
      top_buff[top_index] + (c->hydro.parts - c->top->hydro.parts);
  struct cell_buff *sbuff =
      top_sbuff[top_index] + (c->stars.parts - c->top->stars.parts);
  struct cell_buff *bbuff = top_bbuff[top_index] +
                            (c->black_holes.parts - c->top->black_holes.parts);
  struct cell_buff *gbuff =
      top_gbuff[top_index] + (c->grav.parts - c->top->grav.parts);
  struct cell_buff *sink_buff =
      top_sink_buff[top_index] + (c->sinks.parts - c->top->sinks.parts);

  cell_split(c, c->hydro.parts - s->parts, c->stars.parts - s->sparts,
             c->black_holes.parts - s->bparts,
             c->sinks.parts - s->sinks, buff, sbuff, bbuff, gbuff,
             sink_buff);

  for (int k = 0; k < 8; k++) {
    struct cell *cp = c->progeny[k];

    /* Remove any progeny with zero particles. */
    if (cp->hydro.count == 0 && cp->grav.count == 0 && cp->stars.count == 0 &&
        cp->black_holes.count == 0 && cp->sinks.count == 0) {
      space_recycle(s, cp);
      c->progeny[k] = NULL;
    } else if (depth_remaining > 0) {
      /* Depth budget remains: continue DFS into this child's subtree. */
      space_split_process_subtree(s, top_buff, top_sbuff, top_bbuff,
                                   top_gbuff, top_sink_buff, next, cp, tpid,
                                   depth_remaining - 1);
    } else {
      /* Depth budget exhausted: enqueue for the next BFS level. */
      space_split_frontier_append(next, cp);
    }
  }
}

/**
 * @brief Recursively aggregate per-cell summary data from the leaves up.
 *
 * After the BFS split phase, each leaf cell already has its summary fields
 * (@c h_max, @c ti_*_end_min, multipole, @c maxdepth, etc.) finalised by
 * #space_split_finalize_leaf. This function walks any non-leaf cell's
 * progeny depth-first and reduces their summaries into @c c.
 *
 * Aggregation is the same set of operations that the split branch of the
 * old @c space_split_recursive performed after recursing into its
 * progeny, with one structural difference: the leaf branch is empty here,
 * because leaves have already been finalised.
 *
 * @param s The #space the cell lives in.
 * @param c The cell to aggregate.
 */
static void space_split_aggregate_recursive(struct space *s, struct cell *c) {

  /* Leaves were already finalised during the BFS split phase. */
  if (!c->split) return;

  /* Aggregate progeny first. */
  for (int k = 0; k < 8; k++) {
    if (c->progeny[k] != NULL) space_split_aggregate_recursive(s, c->progeny[k]);
  }

  /* Reduce per-family h_max, time-bin bounds and SF logger across
   * surviving progeny. */
  float h_max = 0.f, h_max_active = 0.f;
  float stars_h_max = 0.f, stars_h_max_active = 0.f;
  float black_holes_h_max = 0.f, black_holes_h_max_active = 0.f;
  float sinks_h_max = 0.f, sinks_h_max_active = 0.f;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_beg_max = 0;
  integertime_t ti_rt_end_min = max_nr_timesteps, ti_rt_beg_max = 0;
  integertime_t ti_rt_min_step_size = max_nr_timesteps;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps, ti_stars_beg_max = 0;
  integertime_t ti_sinks_end_min = max_nr_timesteps, ti_sinks_beg_max = 0;
  integertime_t ti_black_holes_end_min = max_nr_timesteps,
                ti_black_holes_beg_max = 0;
  int maxdepth = 0;

  for (int k = 0; k < 8; k++) {
    const struct cell *cp = c->progeny[k];
    if (cp == NULL) continue;

    h_max = max(h_max, cp->hydro.h_max);
    h_max_active = max(h_max_active, cp->hydro.h_max_active);
    stars_h_max = max(stars_h_max, cp->stars.h_max);
    stars_h_max_active = max(stars_h_max_active, cp->stars.h_max_active);
    black_holes_h_max = max(black_holes_h_max, cp->black_holes.h_max);
    black_holes_h_max_active =
        max(black_holes_h_max_active, cp->black_holes.h_max_active);
    sinks_h_max = max(sinks_h_max, cp->sinks.h_max);
    sinks_h_max_active = max(sinks_h_max_active, cp->sinks.h_max_active);

    ti_hydro_end_min = min(ti_hydro_end_min, cp->hydro.ti_end_min);
    ti_hydro_beg_max = max(ti_hydro_beg_max, cp->hydro.ti_beg_max);
    ti_rt_end_min = min(ti_rt_end_min, cp->rt.ti_rt_end_min);
    ti_rt_beg_max = max(ti_rt_beg_max, cp->rt.ti_rt_beg_max);
    ti_rt_min_step_size =
        min(ti_rt_min_step_size, cp->rt.ti_rt_min_step_size);
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

    maxdepth = max(maxdepth, cp->maxdepth);
  }

  /* Build this cell's multipole from its children (M2M). */
  if (s->with_self_gravity) {

    gravity_reset(c->grav.multipole);

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

    /* Final operation on the CoM and bulk velocity. */
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
    c->grav.multipole->dx_max[0] = 0.f;
    c->grav.multipole->dx_max[1] = 0.f;
    c->grav.multipole->dx_max[2] = 0.f;

    /* Compute the multipole power */
    gravity_multipole_compute_power(&c->grav.multipole->m_pole);
  }

  /* Write the reduced values back to this cell. */
  c->hydro.h_max = h_max;
  c->hydro.h_max_active = h_max_active;
  c->hydro.ti_end_min = ti_hydro_end_min;
  c->hydro.ti_beg_max = ti_hydro_beg_max;
  c->rt.ti_rt_end_min = ti_rt_end_min;
  c->rt.ti_rt_beg_max = ti_rt_beg_max;
  c->rt.ti_rt_min_step_size = ti_rt_min_step_size;
  c->grav.ti_end_min = ti_gravity_end_min;
  c->grav.ti_beg_max = ti_gravity_beg_max;
  c->stars.ti_end_min = ti_stars_end_min;
  c->stars.ti_beg_max = ti_stars_beg_max;
  c->stars.h_max = stars_h_max;
  c->stars.h_max_active = stars_h_max_active;
  c->sinks.ti_end_min = ti_sinks_end_min;
  c->sinks.ti_beg_max = ti_sinks_beg_max;
  c->sinks.h_max = sinks_h_max;
  c->sinks.h_max_active = sinks_h_max_active;
  c->black_holes.ti_end_min = ti_black_holes_end_min;
  c->black_holes.ti_beg_max = ti_black_holes_beg_max;
  c->black_holes.h_max = black_holes_h_max;
  c->black_holes.h_max_active = black_holes_h_max_active;
  c->maxdepth = maxdepth;
  c->owner = -1;
}

/**
 * @brief Decide whether a cell should be split into 8 progeny.
 *
 * The criterion is identical to the one historically used in
 * #space_split_recursive: when self-gravity is on, a cell is split if it
 * contains more than @c space_splitsize gparts; otherwise, it is split when
 * either the hydro or stars particle count exceeds @c space_splitsize.
 *
 * @param s The #space the cell lives in.
 * @param c The #cell to test.
 * @return 1 if @c c should be split, 0 otherwise.
 */
__attribute__((always_inline)) static INLINE int space_split_should_split(
    const struct space *s, const struct cell *c) {

  const int with_self_gravity = s->with_self_gravity;
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;

  if ((with_self_gravity && gcount > space_splitsize) ||
      (!with_self_gravity &&
       (count > space_splitsize || scount > space_splitsize)))
    return 1;

  return 0;
}

/**
 * @brief Allocate and initialise the 8 progeny of a cell that is being split.
 *
 * This sets the parent's @c split flag, draws 8 fresh #cell objects from the
 * sub-cell pool of the executing thread (via #space_getcells) and initialises
 * each of them to be a fresh, empty child cell consistent with the geometry
 * and bookkeeping of @c c. It is a direct extraction of the per-progeny
 * initialisation block previously inlined inside #space_split_recursive.
 *
 * @param s The #space the cell lives in.
 * @param c The parent #cell being split.
 * @param tpid The threadpool ID of the worker that is splitting @c c.
 *        The freshly-allocated children inherit this @c tpid, i.e. cell-pool
 *        ownership follows the worker that performed the split rather than the
 *        original top-level cell's worker.
 */
static void space_split_init_progeny(struct space *s, struct cell *c,
                                     const short int tpid) {

  /* No longer just a leaf. */
  c->split = 1;

  /* Create the cell's progeny. */
  space_getcells(s, 8, c->progeny, tpid);
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
    cp->h_min_allowed = cp->dmin * 0.5 * (1. / kernel_gamma);
    cp->h_max_allowed = cp->dmin * (1. / kernel_gamma);
    if (k & 4) cp->loc[0] += cp->width[0];
    if (k & 2) cp->loc[1] += cp->width[1];
    if (k & 1) cp->loc[2] += cp->width[2];
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
    cp->sinks.h_max = 0.f;
    cp->sinks.h_max_active = 0.f;
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
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
    cell_assign_cell_index(cp, c);
#endif
  }
}

/**
 * @brief Finalise a cell that has been determined to remain a leaf.
 *
 * This collects all the per-particle quantities that the leaf cell exposes to
 * the rest of the code: max smoothing lengths (active and total) for each
 * particle family, time-step bounds, RT time-step bounds, star formation
 * logger contributions, particle-side depth metadata and, when running with
 * self-gravity, the leaf-level multipole built directly from the gparts via
 * #gravity_P2M.
 *
 * It also resets per-particle drift bookkeeping (@c x_diff) and the cell's
 * progeny pointers, sets @c c->maxdepth to @c c->depth and marks the cell as
 * having no current owner. This is a direct extraction of the leaf branch
 * that previously lived inside #space_split_recursive.
 *
 * @param s The #space the cell lives in.
 * @param c The leaf #cell to finalise.
 */
static void space_split_finalize_leaf(struct space *s, struct cell *c) {

  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int bcount = c->black_holes.count;
  const int sink_count = c->sinks.count;
  struct part *parts = c->hydro.parts;
  struct gpart *gparts = c->grav.parts;
  struct spart *sparts = c->stars.parts;
  struct bpart *bparts = c->black_holes.parts;
  struct xpart *xparts = c->hydro.xparts;
  struct sink *sinks = c->sinks.parts;
  struct engine *e = s->e;
  const integertime_t ti_current = e->ti_current;
  const int with_rt = e->policy & engine_policy_rt;

  float h_max = 0.f;
  float h_max_active = 0.f;
  float stars_h_max = 0.f;
  float stars_h_max_active = 0.f;
  float black_holes_h_max = 0.f;
  float black_holes_h_max_active = 0.f;
  float sinks_h_max = 0.f;
  float sinks_h_max_active = 0.f;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_beg_max = 0;
  integertime_t ti_rt_end_min = max_nr_timesteps, ti_rt_beg_max = 0;
  integertime_t ti_rt_min_step_size = max_nr_timesteps;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps, ti_stars_beg_max = 0;
  integertime_t ti_sinks_end_min = max_nr_timesteps, ti_sinks_beg_max = 0;
  integertime_t ti_black_holes_end_min = max_nr_timesteps,
                ti_black_holes_beg_max = 0;

  /* Clear the progeny. */
  bzero(c->progeny, sizeof(struct cell *) * 8);
  c->split = 0;

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
    const timebin_t time_bin_rt = parts[k].rt_time_data.time_bin;
    const integertime_t ti_end = get_integer_time_end(ti_current, time_bin);
    const integertime_t ti_beg = get_integer_time_begin(ti_current, time_bin);

    ti_hydro_end_min = min(ti_hydro_end_min, ti_end);
    ti_hydro_beg_max = max(ti_hydro_beg_max, ti_beg);

    if (with_rt) {
      /* Contrary to other physics, RT doesn't have its own particle type.
       * So collect time step data from particles only when we're running
       * with RT. Otherwise, we may find cells which are active or in
       * impossible timezones. Skipping this check results in cells having
       * RT times = max_nr_timesteps or zero, respecively. */
      const integertime_t ti_rt_end =
          get_integer_time_end(ti_current, time_bin_rt);
      const integertime_t ti_rt_beg =
          get_integer_time_begin(ti_current, time_bin_rt);
      const integertime_t ti_rt_step = get_integer_timestep(time_bin_rt);
      ti_rt_end_min = min(ti_rt_end_min, ti_rt_end);
      ti_rt_beg_max = max(ti_rt_beg_max, ti_rt_beg);
      ti_rt_min_step_size = min(ti_rt_min_step_size, ti_rt_step);
    }

    h_max = max(h_max, parts[k].h);

    if (part_is_active(&parts[k], e))
      h_max_active = max(h_max_active, parts[k].h);

    cell_set_part_h_depth(&parts[k], c);

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
    ti_stars_beg_max = max(ti_stars_beg_max, ti_beg);

    stars_h_max = max(stars_h_max, sparts[k].h);

    if (spart_is_active(&sparts[k], e))
      stars_h_max_active = max(stars_h_max_active, sparts[k].h);

    cell_set_spart_h_depth(&sparts[k], c);

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
    ti_sinks_beg_max = max(ti_sinks_beg_max, ti_beg);

    sinks_h_max = max(sinks_h_max, sinks[k].h);

    if (sink_is_active(&sinks[k], e))
      sinks_h_max_active = max(sinks_h_max_active, sinks[k].h);

    cell_set_sink_h_depth(&sinks[k], c);

    /* Collect SFR from the particles after rebuilt */
    star_formation_logger_log_inactive_sink(&sinks[k], &c->stars.sfh);

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
    ti_black_holes_beg_max = max(ti_black_holes_beg_max, ti_beg);

    black_holes_h_max = max(black_holes_h_max, bparts[k].h);

    if (bpart_is_active(&bparts[k], e))
      black_holes_h_max_active = max(black_holes_h_max_active, bparts[k].h);

    cell_set_bpart_h_depth(&bparts[k], c);

    /* Reset x_diff */
    bparts[k].x_diff[0] = 0.f;
    bparts[k].x_diff[1] = 0.f;
    bparts[k].x_diff[2] = 0.f;
  }

  /* Construct the multipole and the centre of mass*/
  if (s->with_self_gravity) {
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
    c->grav.multipole->dx_max[0] = 0.f;
    c->grav.multipole->dx_max[1] = 0.f;
    c->grav.multipole->dx_max[2] = 0.f;
  }

  /* Set the values for this cell. */
  c->hydro.h_max = h_max;
  c->hydro.h_max_active = h_max_active;
  c->hydro.ti_end_min = ti_hydro_end_min;
  c->hydro.ti_beg_max = ti_hydro_beg_max;
  c->rt.ti_rt_end_min = ti_rt_end_min;
  c->rt.ti_rt_beg_max = ti_rt_beg_max;
  c->rt.ti_rt_min_step_size = ti_rt_min_step_size;
  c->grav.ti_end_min = ti_gravity_end_min;
  c->grav.ti_beg_max = ti_gravity_beg_max;
  c->stars.ti_end_min = ti_stars_end_min;
  c->stars.ti_beg_max = ti_stars_beg_max;
  c->stars.h_max = stars_h_max;
  c->stars.h_max_active = stars_h_max_active;
  c->sinks.ti_end_min = ti_sinks_end_min;
  c->sinks.ti_beg_max = ti_sinks_beg_max;
  c->sinks.h_max = sinks_h_max;
  c->sinks.h_max_active = sinks_h_max_active;
  c->black_holes.ti_end_min = ti_black_holes_end_min;
  c->black_holes.ti_beg_max = ti_black_holes_beg_max;
  c->black_holes.h_max = black_holes_h_max;
  c->black_holes.h_max_active = black_holes_h_max_active;
  c->maxdepth = c->depth;

  /* No runner owns this cell yet. We assign those during scheduling. */
  c->owner = -1;
}

/**
 * @brief Fill the transient top-level split buffers for one top cell.
 *
 * The positions are copied exactly once, the first time a top cell actually
 * decides to split. Descendants later recover their active slice by pointer
 * arithmetic relative
 * to @c c->top and reuse the already-permuted buffers just like the legacy DFS
 * implementation did when it passed slices down recursively.
 *
 * @param c The top-level #cell owning the particle ranges.
 * @param buff Hydro particle buffer.
 * @param sbuff Stars particle buffer.
 * @param bbuff Black-hole particle buffer.
 * @param gbuff Gravity particle buffer.
 * @param sink_buff Sink particle buffer.
 */
static void space_split_fill_top_buffers(const struct cell *c,
                                         struct cell_buff *buff,
                                         struct cell_buff *sbuff,
                                         struct cell_buff *bbuff,
                                         struct cell_buff *gbuff,
                                         struct cell_buff *sink_buff) {

  for (int k = 0; k < c->hydro.count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->hydro.parts[k].time_bin == time_bin_inhibited)
      error("Inhibited particle present in space_split()");
    if (c->hydro.parts[k].time_bin == time_bin_not_created)
      error("Extra particle present in space_split()");
#endif
    buff[k].x[0] = c->hydro.parts[k].x[0];
    buff[k].x[1] = c->hydro.parts[k].x[1];
    buff[k].x[2] = c->hydro.parts[k].x[2];
  }

  for (int k = 0; k < c->grav.count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->grav.parts[k].time_bin == time_bin_inhibited)
      error("Inhibited particle present in space_split()");
    if (c->grav.parts[k].time_bin == time_bin_not_created)
      error("Extra particle present in space_split()");
#endif
    gbuff[k].x[0] = c->grav.parts[k].x[0];
    gbuff[k].x[1] = c->grav.parts[k].x[1];
    gbuff[k].x[2] = c->grav.parts[k].x[2];
  }

  for (int k = 0; k < c->stars.count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->stars.parts[k].time_bin == time_bin_inhibited)
      error("Inhibited particle present in space_split()");
    if (c->stars.parts[k].time_bin == time_bin_not_created)
      error("Extra particle present in space_split()");
#endif
    sbuff[k].x[0] = c->stars.parts[k].x[0];
    sbuff[k].x[1] = c->stars.parts[k].x[1];
    sbuff[k].x[2] = c->stars.parts[k].x[2];
  }

  for (int k = 0; k < c->black_holes.count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->black_holes.parts[k].time_bin == time_bin_inhibited)
      error("Inhibited particle present in space_split()");
    if (c->black_holes.parts[k].time_bin == time_bin_not_created)
      error("Extra particle present in space_split()");
#endif
    bbuff[k].x[0] = c->black_holes.parts[k].x[0];
    bbuff[k].x[1] = c->black_holes.parts[k].x[1];
    bbuff[k].x[2] = c->black_holes.parts[k].x[2];
  }

  for (int k = 0; k < c->sinks.count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (c->sinks.parts[k].time_bin == time_bin_inhibited)
      error("Inhibited particle present in space_split()");
    if (c->sinks.parts[k].time_bin == time_bin_not_created)
      error("Extra particle present in space_split()");
#endif
    sink_buff[k].x[0] = c->sinks.parts[k].x[0];
    sink_buff[k].x[1] = c->sinks.parts[k].x[1];
    sink_buff[k].x[2] = c->sinks.parts[k].x[2];
  }
}


/**
 * @brief #threadpool mapper function for one BFS level.
 *
 * Each chunk of work is a contiguous range of cells from the
 * @c current frontier. Each cell processes a limited number of DFS levels
 * via #space_split_process_subtree before surviving children are appended
 * to the @c next frontier.
 *
 * @param map_data Pointer to the start of a chunk of #cell pointers.
 * @param num_cells Number of cells in this chunk.
 * @param extra_data Pointer to an array containing the #space, the transient
 *        top-level split buffers and the next-level frontier.
 */
static void space_split_frontier_mapper(void *map_data, int num_cells,
                                        void *extra_data) {

  void **mapper_data = (void **)extra_data;
  struct space *s = (struct space *)mapper_data[0];
  struct cell_buff **top_buff = (struct cell_buff **)mapper_data[1];
  struct cell_buff **top_sbuff = (struct cell_buff **)mapper_data[2];
  struct cell_buff **top_bbuff = (struct cell_buff **)mapper_data[3];
  struct cell_buff **top_gbuff = (struct cell_buff **)mapper_data[4];
  struct cell_buff **top_sink_buff = (struct cell_buff **)mapper_data[5];
  struct space_split_frontier *next =
      (struct space_split_frontier *)mapper_data[6];
  struct cell **cells = (struct cell **)map_data;
  const short int tpid = threadpool_gettid();

  for (int i = 0; i < num_cells; i++)
    space_split_process_subtree(s, top_buff, top_sbuff, top_bbuff,
                                top_gbuff, top_sink_buff, next, cells[i],
                                tpid, dfs_levels_at_a_time - 1);
}



/**
 * @brief #threadpool mapper function for the upward aggregation pass.
 *
 * Each chunk is a contiguous range of top-level cell indices. The
 * aggregation walks each top cell's subtree depth-first via
 * #space_split_aggregate_recursive and then folds in the global
 * top-level multipole statistics (@c min_a_grav, @c max_softening,
 * @c max_mpole_power) using atomic min/max reductions.
 *
 * @param map_data Pointer to the start of a chunk of cell indices.
 * @param num_cells Number of indices in this chunk.
 * @param extra_data Pointer to the parent #space.
 */
static void space_split_aggregate_mapper(void *map_data, int num_cells,
                                         void *extra_data) {

  struct space *s = (struct space *)extra_data;
  struct cell *cells_top = s->cells_top;
  int *local_cells_with_particles = (int *)map_data;

  float min_a_grav = FLT_MAX;
  float max_softening = 0.f;
  float max_mpole_power[SELF_GRAVITY_MULTIPOLE_ORDER + 1] = {0.f};

  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells_with_particles[ind]];

    space_split_aggregate_recursive(s, c);

    if (s->with_self_gravity) {
      min_a_grav =
          min(min_a_grav, c->grav.multipole->m_pole.min_old_a_grav_norm);
      max_softening =
          max(max_softening, c->grav.multipole->m_pole.max_softening);

      for (int n = 0; n < SELF_GRAVITY_MULTIPOLE_ORDER + 1; ++n)
        max_mpole_power[n] =
            max(max_mpole_power[n], c->grav.multipole->m_pole.power[n]);
    }

    /* The deepest cell in this top tree determines the contribution to
     * the global max depth. */
    atomic_max(&s->maxdepth, c->maxdepth);
  }

#ifdef SWIFT_DEBUG_CHECKS
  for (int ind = 0; ind < num_cells; ind++) {
    int depth = 0;
    const struct cell *c = &cells_top[local_cells_with_particles[ind]];
    if (!checkCellhdxmax(c, &depth)) message("    at cell depth %d", depth);
  }
#endif

  atomic_min_f(&s->min_a_grav, min_a_grav);
  atomic_max_f(&s->max_softening, max_softening);
  for (int n = 0; n < SELF_GRAVITY_MULTIPOLE_ORDER + 1; ++n)
    atomic_max_f(&s->max_mpole_power[n], max_mpole_power[n]);
}

/**
 * @brief Split particles between cells of a hierarchy.
 *
 * BFS-parallel implementation with limited DFS per frontier entry:
 *
 * 1. <em>BFS split phase.</em> A frontier of cells (initially the local
 *    top-level cells with particles) is processed level by level using
 *    #threadpool_map. At each level each cell processes up to
 *    @c dfs_levels_at_a_time levels depth-first via
 *    #space_split_process_subtree before survivors re-enter the frontier.
 *    This combines the cache-locality benefits of DFS within a small subtree
 *    with the work-stealing benefits of BFS at the frontier: a deep tree
 *    distributing its cells across many BFS rounds is broken into batches of
 *    @c dfs_levels_at_a_time levels, keeping a thread's data hot for that
 *    batch and then re-balancing at the frontier.
 *
 * 2. <em>Aggregation phase.</em> Once all splitting is done, each top cell's
 *    subtree is walked depth-first to combine per-family @c h_max, time-bin
 *    bounds and multipoles from the leaves up to the top via
 *    #space_split_aggregate_recursive. This phase is parallel over top cells;
 *    each subtree is independent.
 *
 * Each top cell that actually splits owns transient split buffers that are
 * allocated and filled exactly once, on demand. Descendants recover their
 * current buffer slice from @c c->top and the particle-array offsets within
 * that top cell, matching the legacy DFS buffer semantics without carrying
 * extra data in the frontier. The BFS
 * frontier arrays are allocated once and sized to the worst case
 * (@c 8 * size_of_current_frontier) before each level, so
 * #space_split_frontier_append needs only an atomic counter increment to
 * reserve a slot.
 *
 * @param s The #space.
 * @param verbose Are we talkative ?
 */
void space_split(struct space *s, int verbose) {

  const ticks tic = getticks();

  s->min_a_grav = FLT_MAX;
  s->max_softening = 0.f;
  bzero(s->max_mpole_power, (SELF_GRAVITY_MULTIPOLE_ORDER + 1) * sizeof(float));

  const int nr_local_cells = s->nr_local_cells_with_particles;

  if (nr_local_cells > 0) {

    struct threadpool *tp = &s->e->threadpool;
    const int nr_top_cells = s->nr_cells;

    struct cell_buff **top_buff = (struct cell_buff **)swift_malloc(
        "top_split_buff", sizeof(struct cell_buff *) * nr_top_cells);
    struct cell_buff **top_sbuff = (struct cell_buff **)swift_malloc(
        "top_split_sbuff", sizeof(struct cell_buff *) * nr_top_cells);
    struct cell_buff **top_bbuff = (struct cell_buff **)swift_malloc(
        "top_split_bbuff", sizeof(struct cell_buff *) * nr_top_cells);
    struct cell_buff **top_gbuff = (struct cell_buff **)swift_malloc(
        "top_split_gbuff", sizeof(struct cell_buff *) * nr_top_cells);
    struct cell_buff **top_sink_buff = (struct cell_buff **)swift_malloc(
        "top_split_sink_buff", sizeof(struct cell_buff *) * nr_top_cells);

    if (top_buff == NULL || top_sbuff == NULL || top_bbuff == NULL ||
        top_gbuff == NULL || top_sink_buff == NULL)
      error("Failed to allocate top-level split buffer pointer arrays.");

    bzero(top_buff, sizeof(struct cell_buff *) * nr_top_cells);
    bzero(top_sbuff, sizeof(struct cell_buff *) * nr_top_cells);
    bzero(top_bbuff, sizeof(struct cell_buff *) * nr_top_cells);
    bzero(top_gbuff, sizeof(struct cell_buff *) * nr_top_cells);
    bzero(top_sink_buff, sizeof(struct cell_buff *) * nr_top_cells);

    /* Two frontier buffers, swapped each level. One holds the cells of the
     * current BFS level while the other accumulates their surviving children
     * for the next level. Seed the current-level frontier with the local
     * top-level cells with particles. */
    struct space_split_frontier this_level_frontier = {NULL, 0, 0};
    struct space_split_frontier next_level_frontier = {NULL, 0, 0};
    struct space_split_frontier *this_level = &this_level_frontier;
    struct space_split_frontier *next_level = &next_level_frontier;

    space_split_frontier_ensure_capacity(this_level, nr_local_cells);
    for (int ind = 0; ind < nr_local_cells; ind++) {
      struct cell *c =
          &s->cells_top[s->local_cells_with_particles_top[ind]];
      this_level->cells[ind] = c;
    }
    this_level->count = nr_local_cells;

    void *frontier_mapper_data[7] = {s,
                                     top_buff,
                                     top_sbuff,
                                     top_bbuff,
                                     top_gbuff,
                                     top_sink_buff,
                                     next_level};

    /* Drive the BFS until the frontier is empty. Each cell processes up to
     * @c dfs_levels_at_a_time levels depth-first before survivors re-enter
     * the frontier at the next level. */
    const ticks bfs_tic = getticks();
    while (this_level->count > 0) {

      /* Pre-size the next frontier to the worst case (each cell can
       * spawn up to 8 surviving children after its DFS batch). */
      space_split_frontier_ensure_capacity(next_level, 8 * this_level->count);
      next_level->count = 0;

      frontier_mapper_data[6] = next_level;

      threadpool_map(tp, space_split_frontier_mapper, this_level->cells,
                     this_level->count, sizeof(struct cell *),
                     threadpool_auto_chunk_size, frontier_mapper_data);

      /* Swap frontiers. */
      struct space_split_frontier *tmp = this_level;
      this_level = next_level;
      next_level = tmp;
    }
    const ticks bfs_toc = getticks();

    /* Release the transient top-level split buffers and the frontier
     * storage. */
    for (int ind = 0; ind < nr_local_cells; ind++) {
      const int top_index = s->local_cells_with_particles_top[ind];
      if (top_buff[top_index] != NULL)
        swift_free("top_split_buff", top_buff[top_index]);
      if (top_sbuff[top_index] != NULL)
        swift_free("top_split_sbuff", top_sbuff[top_index]);
      if (top_bbuff[top_index] != NULL)
        swift_free("top_split_bbuff", top_bbuff[top_index]);
      if (top_gbuff[top_index] != NULL)
        swift_free("top_split_gbuff", top_gbuff[top_index]);
      if (top_sink_buff[top_index] != NULL)
        swift_free("top_split_sink_buff", top_sink_buff[top_index]);
    }

    swift_free("top_split_buff", top_buff);
    swift_free("top_split_sbuff", top_sbuff);
    swift_free("top_split_bbuff", top_bbuff);
    swift_free("top_split_gbuff", top_gbuff);
    swift_free("top_split_sink_buff", top_sink_buff);

    if (this_level_frontier.cells != NULL)
      swift_free("split_frontier", this_level_frontier.cells);
    if (next_level_frontier.cells != NULL)
      swift_free("split_frontier", next_level_frontier.cells);

    /* Aggregation phase: walk each top tree from leaves up to combine
     * per-family h_max, time-bin bounds and multipoles. */
    const ticks aggregate_tic = getticks();
    threadpool_map(tp, space_split_aggregate_mapper,
                   s->local_cells_with_particles_top, nr_local_cells,
                   sizeof(int), threadpool_auto_chunk_size, s);
    const ticks aggregate_toc = getticks();

    if (verbose) {
      message("BFS split loop: %.3f %s.", clocks_from_ticks(bfs_toc - bfs_tic),
              clocks_getunit());
      message("Aggregate pass: %.3f %s.",
              clocks_from_ticks(aggregate_toc - aggregate_tic),
              clocks_getunit());
    }
  }

  if (verbose) {
    message("Max tree depth after split: %d", s->maxdepth);
    message("Have %d cells including subcells (cell footprint: %zd MB)",
            s->tot_cells, s->tot_cells * sizeof(struct cell) / (1024 * 1024));
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
  }
}
