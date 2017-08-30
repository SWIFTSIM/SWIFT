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
#ifndef SWIFT_CELL_H
#define SWIFT_CELL_H

/* Config parameters. */
#include "../config.h"

/* Includes. */
#include <stddef.h>

/* Local includes. */
#include "align.h"
#include "kernel_hydro.h"
#include "lock.h"
#include "multipole.h"
#include "part.h"
#include "space.h"
#include "task.h"
#include "timeline.h"

/* Avoid cyclic inclusions */
struct engine;
struct scheduler;

/* Max tag size set to 2^29 to take into account some MPI implementations
 * that use 2^31 as the upper bound on MPI tags and the fact that
 * cell_next_tag is multiplied by 2 when passed to an MPI function.
 * The maximum was lowered by a further factor of 2 to be on the safe side.*/
#define cell_max_tag (1 << 29)

#define cell_align 128

/* Global variables. */
extern int cell_next_tag;

/* Struct to temporarily buffer the particle locations and bin id. */
struct cell_buff {
  double x[3];
  int ind;
} SWIFT_STRUCT_ALIGN;

/* Mini struct to link cells to tasks. Used as a linked list. */
struct link {

  /* The task pointer. */
  struct task *t;

  /* The next pointer. */
  struct link *next;
};

/* Packed cell. */
struct pcell {

  /* Stats on this cell's particles. */
  double h_max;
  integertime_t ti_end_min, ti_end_max, ti_beg_max, ti_old_part, ti_old_gpart;

  /* Number of particles in this cell. */
  int count, gcount, scount;

  /* tag used for MPI communication. */
  int tag;

  /* Relative indices of the cell's progeny. */
  int progeny[8];

} SWIFT_STRUCT_ALIGN;

/**
 * @brief Cell within the tree structure.
 *
 * Contains particles, links to tasks, a multipole object and counters.
 */
struct cell {

  /*! The cell location on the grid. */
  double loc[3];

  /*! The cell dimensions. */
  double width[3];

  /*! Max smoothing length in this cell. */
  double h_max;

  /*! This cell's multipole. */
  struct gravity_tensors *multipole;

  /*! Linking pointer for "memory management". */
  struct cell *next;

  /*! Pointer to the #part data. */
  struct part *parts;

  /*! Pointer to the #xpart data. */
  struct xpart *xparts;

  /*! Pointer to the #gpart data. */
  struct gpart *gparts;

  /*! Pointer to the #spart data. */
  struct spart *sparts;

  /*! Pointer for the sorted indices. */
  struct entry *sort[13];

  /*! Pointers to the next level of cells. */
  struct cell *progeny[8];

  /*! Parent cell. */
  struct cell *parent;

  /*! Super cell, i.e. the highest-level parent cell that has pair/self tasks */
  struct cell *super;

  /*! The task computing this cell's sorts. */
  struct task *sorts;

  /*! Linked list of the tasks computing this cell's hydro density. */
  struct link *density;

  /* Linked list of the tasks computing this cell's hydro gradients. */
  struct link *gradient;

  /*! Linked list of the tasks computing this cell's hydro forces. */
  struct link *force;

  /*! Linked list of the tasks computing this cell's gravity forces. */
  struct link *grav;

  /*! The multipole initialistation task */
  struct task *init_grav;

  /*! Dependency implicit task for the ghost  (in->ghost->out)*/
  struct task *ghost_in;

  /*! Dependency implicit task for the ghost  (in->ghost->out)*/
  struct task *ghost_out;

  /*! The ghost task itself */
  struct task *ghost;

  /*! The extra ghost task for complex hydro schemes */
  struct task *extra_ghost;

  /*! The drift task for parts */
  struct task *drift_part;

  /*! The drift task for gparts */
  struct task *drift_gpart;

  /*! The first kick task */
  struct task *kick1;

  /*! The second kick task */
  struct task *kick2;

  /*! The task to compute time-steps */
  struct task *timestep;

  /*! Task linking the FFT mesh to the rest of gravity tasks */
  struct task *grav_ghost[2];

  /*! Task computing long range non-periodic gravity interactions */
  struct task *grav_long_range;

  /*! Task propagating the multipole to the particles */
  struct task *grav_down;

  /*! Task for cooling */
  struct task *cooling;

  /*! Task for source terms */
  struct task *sourceterms;

#ifdef WITH_MPI

  /* Task receiving data (positions). */
  struct task *recv_xv;

  /* Task receiving data (density). */
  struct task *recv_rho;

  /* Task receiving data (gradient). */
  struct task *recv_gradient;

  /* Task receiving data (time-step). */
  struct task *recv_ti;

  /* Linked list for sending data (positions). */
  struct link *send_xv;

  /* Linked list for sending data (density). */
  struct link *send_rho;

  /* Linked list for sending data (gradient). */
  struct link *send_gradient;

  /* Linked list for sending data (time-step). */
  struct link *send_ti;

  /*! Bit mask of the proxies this cell is registered with. */
  unsigned long long int sendto;

  /*! Pointer to this cell's packed representation. */
  struct pcell *pcell;

  /*! Size of the packed representation */
  int pcell_size;

  /*! MPI tag associated with this cell */
  int tag;

#endif

  /*! Minimum end of (integer) time step in this cell. */
  integertime_t ti_end_min;

  /*! Maximum end of (integer) time step in this cell. */
  integertime_t ti_end_max;

  /*! Maximum beginning of (integer) time step in this cell. */
  integertime_t ti_beg_max;

  /*! Last (integer) time the cell's part were drifted forward in time. */
  integertime_t ti_old_part;

  /*! Last (integer) time the cell's gpart were drifted forward in time. */
  integertime_t ti_old_gpart;

  /*! Last (integer) time the cell's multipole was drifted forward in time. */
  integertime_t ti_old_multipole;

  /*! Minimum dimension, i.e. smallest edge of this cell (min(width)). */
  float dmin;

  /*! Maximum particle movement in this cell since the last sort. */
  float dx_max_sort;

  /*! Maximum part movement in this cell since last construction. */
  float dx_max_part;

  /*! Maximum gpart movement in this cell since last construction. */
  float dx_max_gpart;

  /*! Nr of #part in this cell. */
  int count;

  /*! Nr of #gpart in this cell. */
  int gcount;

  /*! Nr of #spart in this cell. */
  int scount;

  /*! Bit-mask indicating the sorted directions */
  unsigned int sorted;

  /*! Spin lock for various uses (#part case). */
  swift_lock_type lock;

  /*! Spin lock for various uses (#gpart case). */
  swift_lock_type glock;

  /*! Spin lock for various uses (#multipole case). */
  swift_lock_type mlock;

  /*! Spin lock for various uses (#spart case). */
  swift_lock_type slock;

  /*! ID of the previous owner, e.g. runner. */
  int owner;

  /*! Number of #part updated in this cell. */
  int updated;

  /*! Number of #gpart updated in this cell. */
  int g_updated;

  /*! Number of #spart updated in this cell. */
  int s_updated;

  /*! ID of the node this cell lives on. */
  int nodeID;

  /*! Is the #part data of this cell being used in a sub-cell? */
  int hold;

  /*! Is the #gpart data of this cell being used in a sub-cell? */
  int ghold;

  /*! Is the #multipole data of this cell being used in a sub-cell? */
  int mhold;

  /*! Is the #spart data of this cell being used in a sub-cell? */
  int shold;

  /*! Values of dx_max before the drifts, used for sub-cell tasks. */
  float dx_max_old;

  /*! Values of h_max before the drifts, used for sub-cell tasks. */
  float h_max_old;

  /*! Values of dx_max_sort before the drifts, used for sub-cell tasks. */
  float dx_max_sort_old;

  /*! Bit mask of sort directions that will be needed in the next timestep. */
  unsigned int requires_sorts;

  /*! Bit mask of sorts that need to be computed for this cell. */
  unsigned int do_sort;

  /*! Number of tasks that are associated with this cell. */
  short int nr_tasks;

  /*! The depth of this cell in the tree. */
  char depth;

  /*! Is this cell split ? */
  char split;

  /*! The maximal depth of this cell and its progenies */
  char maxdepth;

  /*! Does this cell need to be drifted (hydro)? */
  char do_drift;

  /*! Do any of this cell's sub-cells need to be drifted (hydro)? */
  char do_sub_drift;

  /*! Does this cell need to be drifted (gravity)? */
  char do_grav_drift;

  /*! Do any of this cell's sub-cells need to be drifted (gravity)? */
  char do_grav_sub_drift;

  /*! Do any of this cell's sub-cells need to be sorted? */
  char do_sub_sort;

#ifdef SWIFT_DEBUG_CHECKS
  /*! Last (integer) time the cell's sort arrays were updated. */
  integertime_t ti_sort;

  /*! The list of tasks that have been executed on this cell */
  char tasks_executed[64];

  /*! The list of sub-tasks that have been executed on this cell */
  char subtasks_executed[64];
#endif

} SWIFT_STRUCT_ALIGN;

/* Convert cell location to ID. */
#define cell_getid(cdim, i, j, k) \
  ((int)(k) + (cdim)[2] * ((int)(j) + (cdim)[1] * (int)(i)))

/* Function prototypes. */
void cell_split(struct cell *c, ptrdiff_t parts_offset, ptrdiff_t sparts_offset,
                struct cell_buff *buff, struct cell_buff *sbuff,
                struct cell_buff *gbuff);
void cell_sanitize(struct cell *c, int treated);
int cell_locktree(struct cell *c);
void cell_unlocktree(struct cell *c);
int cell_glocktree(struct cell *c);
void cell_gunlocktree(struct cell *c);
int cell_mlocktree(struct cell *c);
void cell_munlocktree(struct cell *c);
int cell_slocktree(struct cell *c);
void cell_sunlocktree(struct cell *c);
int cell_pack(struct cell *c, struct pcell *pc);
int cell_unpack(struct pcell *pc, struct cell *c, struct space *s);
int cell_pack_ti_ends(struct cell *c, integertime_t *ti_ends);
int cell_unpack_ti_ends(struct cell *c, integertime_t *ti_ends);
int cell_getsize(struct cell *c);
int cell_link_parts(struct cell *c, struct part *parts);
int cell_link_gparts(struct cell *c, struct gpart *gparts);
int cell_link_sparts(struct cell *c, struct spart *sparts);
void cell_convert_hydro(struct cell *c, void *data);
void cell_clean_links(struct cell *c, void *data);
void cell_make_multipoles(struct cell *c, integertime_t ti_current);
void cell_check_multipole(struct cell *c, void *data);
void cell_clean(struct cell *c);
void cell_check_part_drift_point(struct cell *c, void *data);
void cell_check_gpart_drift_point(struct cell *c, void *data);
void cell_check_multipole_drift_point(struct cell *c, void *data);
void cell_reset_task_counters(struct cell *c);
int cell_unskip_tasks(struct cell *c, struct scheduler *s);
void cell_set_super(struct cell *c, struct cell *super);
void cell_drift_part(struct cell *c, const struct engine *e, int force);
void cell_drift_gpart(struct cell *c, const struct engine *e, int force);
void cell_drift_multipole(struct cell *c, const struct engine *e);
void cell_drift_all_multipoles(struct cell *c, const struct engine *e);
void cell_check_timesteps(struct cell *c);
void cell_store_pre_drift_values(struct cell *c);
void cell_activate_subcell_tasks(struct cell *ci, struct cell *cj,
                                 struct scheduler *s);
void cell_activate_subcell_grav_tasks(struct cell *ci, struct cell *cj,
                                      struct scheduler *s);
void cell_activate_subcell_external_grav_tasks(struct cell *ci,
                                               struct scheduler *s);
void cell_activate_drift_part(struct cell *c, struct scheduler *s);
void cell_activate_drift_gpart(struct cell *c, struct scheduler *s);
void cell_activate_sorts(struct cell *c, int sid, struct scheduler *s);
void cell_clear_drift_flags(struct cell *c, void *data);
void cell_set_super_mapper(void *map_data, int num_elements, void *extra_data);

/* Inlined functions (for speed). */

/**
 * @brief Can a sub-pair hydro task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_can_recurse_in_pair_task(
    const struct cell *c) {

  /* Is the cell split ? */
  /* If so, is the cut-off radius plus the max distance the parts have moved */
  /* smaller than the sub-cell sizes ? */
  /* Note: We use the _old values as these might have been updated by a drift */
  return c->split &&
         ((kernel_gamma * c->h_max_old + c->dx_max_old) < 0.5f * c->dmin);
}

/**
 * @brief Can a sub-self hydro task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_can_recurse_in_self_task(
    const struct cell *c) {

  /* Is the cell split ? */
  /* Note: No need for more checks here as all the sub-pairs and sub-self */
  /* operations will be executed. So no need for the particle to be at exactly
   */
  /* the right place. */
  return c->split;
}

/**
 * @brief Can a pair task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_can_split_pair_task(
    const struct cell *c) {

  /* Is the cell split ? */
  /* If so, is the cut-off radius with some leeway smaller than */
  /* the sub-cell sizes ? */
  /* Note that since tasks are create after a rebuild no need to take */
  /* into account any part motion (i.e. dx_max == 0 here) */
  return c->split && (space_stretch * kernel_gamma * c->h_max < 0.5f * c->dmin);
}

/**
 * @brief Can a self task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_can_split_self_task(
    const struct cell *c) {

  /* Is the cell split ? */
  /* Note: No need for more checks here as all the sub-pairs and sub-self */
  /* tasks will be created. So no need to check for h_max */
  return c->split && (space_stretch * kernel_gamma * c->h_max < 0.5f * c->dmin);
}

/**
 * @brief Have particles in a pair of cells moved too much and require a rebuild
 * ?
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
__attribute__((always_inline)) INLINE static int cell_need_rebuild_for_pair(
    const struct cell *ci, const struct cell *cj) {

  /* Is the cut-off radius plus the max distance the parts in both cells have */
  /* moved larger than the cell size ? */
  /* Note ci->dmin == cj->dmin */
  return (kernel_gamma * max(ci->h_max, cj->h_max) + ci->dx_max_part +
              cj->dx_max_part >
          cj->dmin);
}

#endif /* SWIFT_CELL_H */
