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

/**
 * @brief Packed cell for information correct at rebuild time.
 *
 * Contains all the information for a tree walk in a non-local cell.
 */
struct pcell {

  /*! Hydro variables */
  struct {

    /*! Maximal smoothing length. */
    double h_max;

    /*! Minimal integer end-of-timestep in this cell for hydro tasks */
    integertime_t ti_end_min;

    /*! Maximal integer end-of-timestep in this cell for hydro tasks */
    integertime_t ti_end_max;

    /*! Maximal integer beginning-of-timestep in this cell for hydro tasks */
    integertime_t ti_beg_max;

    /*! Integer time of the last drift of the #part in this cell */
    integertime_t ti_old_part;

    /*! Number of #part in this cell. */
    int count;

  } hydro;

  /*! Gravity variables */
  struct {

    /*! This cell's gravity-related tensors */
    struct multipole m_pole;

    /*! Centre of mass. */
    double CoM[3];

    /*! Centre of mass at rebuild time. */
    double CoM_rebuild[3];

    /*! Upper limit of the CoM<->gpart distance. */
    double r_max;

    /*! Upper limit of the CoM<->gpart distance at last rebuild. */
    double r_max_rebuild;

    /*! Minimal integer end-of-timestep in this cell for gravity tasks */
    integertime_t ti_end_min;

    /*! Maximal integer end-of-timestep in this cell for gravity tasks */
    integertime_t ti_end_max;

    /*! Maximal integer beginning-of-timestep in this cell for gravity tasks */
    integertime_t ti_beg_max;

    /*! Integer time of the last drift of the #gpart in this cell */
    integertime_t ti_old_part;

    /*! Integer time of the last drift of the #multipole in this cell */
    integertime_t ti_old_multipole;

    /*! Number of #gpart in this cell. */
    int count;

  } grav;

  /*! Stars variables */
  struct {

    /*! Number of #spart in this cell. */
    int count;

  } stars;

  /*! Relative indices of the cell's progeny. */
  int progeny[8];

#ifdef SWIFT_DEBUG_CHECKS
  /* Cell ID (for debugging) */
  int cellID;
#endif

} SWIFT_STRUCT_ALIGN;

/**
 * @brief Cell information at the end of a time-step.
 */
struct pcell_step {

  /*! Hydro variables */
  struct {

    /*! Minimal integer end-of-timestep in this cell (hydro) */
    integertime_t ti_end_min;

    /*! Minimal integer end-of-timestep in this cell (hydro) */
    integertime_t ti_end_max;

    /*! Maximal distance any #part has travelled since last rebuild */
    float dx_max_part;

  } hydro;

  /*! Grav variables */
  struct {

    /*! Minimal integer end-of-timestep in this cell (gravity) */
    integertime_t ti_end_min;

    /*! Minimal integer end-of-timestep in this cell (gravity) */
    integertime_t ti_end_max;

  } grav;

  /*! Stars variables */
  struct {
  } stars;
};

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

  /*! Linking pointer for "memory management". */
  struct cell *next;

  /*! Pointers to the next level of cells. */
  struct cell *progeny[8];

  /*! Parent cell. */
  struct cell *parent;

  /*! Super cell, i.e. the highest-level parent cell with *any* task */
  struct cell *super;

  /*! Hydro variables */
  struct {

    /*! Pointer to the #part data. */
    struct part *parts;

    /*! Pointer to the #xpart data. */
    struct xpart *xparts;

    /*! Pointer for the sorted indices. */
    struct entry *sort[13];

    /*! Super cell, i.e. the highest-level parent cell that has a hydro
     * pair/self tasks */
    struct cell *super;

    /*! Last (integer) time the cell's part were drifted forward in time. */
    integertime_t ti_old_part;

    /*! Maximum part movement in this cell since last construction. */
    float dx_max_part;

    /*! Maximum particle movement in this cell since the last sort. */
    float dx_max_sort;

    /*! Max smoothing length in this cell. */
    double h_max;

    /*! Minimum end of (integer) time step in this cell for hydro tasks. */
    integertime_t ti_end_min;

    /*! Maximum end of (integer) time step in this cell for hydro tasks. */
    integertime_t ti_end_max;

    /*! Maximum beginning of (integer) time step in this cell for hydro tasks.
     */
    integertime_t ti_beg_max;

    /*! Nr of #part in this cell. */
    int count;

    /*! Spin lock for various uses (#part case). */
    swift_lock_type lock;

    /*! Number of #part updated in this cell. */
    int updated;

    /*! Is the #part data of this cell being used in a sub-cell? */
    int hold;

    /*! Values of h_max before the drifts, used for sub-cell tasks. */
    float h_max_old;

    /*! Values of dx_max before the drifts, used for sub-cell tasks. */
    float dx_max_part_old;

    /*! Values of dx_max_sort before the drifts, used for sub-cell tasks. */
    float dx_max_sort_old;

    /*! Bit mask of sort directions that will be needed in the next timestep. */
    unsigned int requires_sorts;

    /*! Bit mask of sorts that need to be computed for this cell. */
    unsigned int do_sort;

    /*! Does this cell need to be drifted (hydro)? */
    char do_drift;

    /*! Do any of this cell's sub-cells need to be drifted (hydro)? */
    char do_sub_drift;

    /*! Do any of this cell's sub-cells need to be sorted? */
    char do_sub_sort;

    /*! Bit-mask indicating the sorted directions */
    unsigned int sorted;

    /*! The task computing this cell's sorts. */
    struct task *sorts;

    /*! The drift task for parts */
    struct task *drift;

    /*! Linked list of the tasks computing this cell's hydro density. */
    struct link *density;

    /* Linked list of the tasks computing this cell's hydro gradients. */
    struct link *gradient;

    /*! Linked list of the tasks computing this cell's hydro forces. */
    struct link *force;

    /*! Dependency implicit task for the ghost  (in->ghost->out)*/
    struct task *ghost_in;

    /*! Dependency implicit task for the ghost  (in->ghost->out)*/
    struct task *ghost_out;

    /*! The ghost task itself */
    struct task *ghost;

    /*! The extra ghost task for complex hydro schemes */
    struct task *extra_ghost;

    /*! The task to end the force calculation */
    struct task *end_force;

    /*! Task for cooling */
    struct task *cooling;

#ifdef SWIFT_DEBUG_CHECKS

    /*! Last (integer) time the cell's sort arrays were updated. */
    integertime_t ti_sort;

#endif

  } hydro;

  /*! Grav variables */
  struct {

    /*! Pointer to the #gpart data. */
    struct gpart *parts;

    /*! This cell's multipole. */
    struct gravity_tensors *multipole;

    /*! Super cell, i.e. the highest-level parent cell that has a grav pair/self
     * tasks */
    struct cell *super;

    /*! Minimum end of (integer) time step in this cell for gravity tasks. */
    integertime_t ti_end_min;

    /*! Maximum end of (integer) time step in this cell for gravity tasks. */
    integertime_t ti_end_max;

    /*! Maximum beginning of (integer) time step in this cell for gravity tasks.
     */
    integertime_t ti_beg_max;

    /*! Last (integer) time the cell's gpart were drifted forward in time. */
    integertime_t ti_old_part;

    /*! Last (integer) time the cell's multipole was drifted forward in time. */
    integertime_t ti_old_multipole;

    /*! Nr of #gpart in this cell. */
    int count;

    /*! Spin lock for various uses (#gpart case). */
    swift_lock_type plock;

    /*! Spin lock for various uses (#multipole case). */
    swift_lock_type mlock;

    /*! Number of #gpart updated in this cell. */
    int updated;

    /*! Is the #gpart data of this cell being used in a sub-cell? */
    int phold;

    /*! Is the #multipole data of this cell being used in a sub-cell? */
    int mhold;

    /*! Does this cell need to be drifted (gravity)? */
    char do_drift;

    /*! Do any of this cell's sub-cells need to be drifted (gravity)? */
    char do_sub_drift;

    /*! The drift task for gparts */
    struct task *drift;

    /*! Linked list of the tasks computing this cell's gravity forces. */
    struct link *grav;

    /*! Linked list of the tasks computing this cell's gravity M-M forces. */
    struct link *mm;

    /*! The multipole initialistation task */
    struct task *init;

    /*! Implicit task for the gravity initialisation */
    struct task *init_out;

    /*! Task computing long range non-periodic gravity interactions */
    struct task *long_range;

    /*! Implicit task for the down propagation */
    struct task *down_in;

    /*! Task propagating the mesh forces to the particles */
    struct task *mesh;

    /*! Task propagating the multipole to the particles */
    struct task *down;

    /*! Number of M-M tasks that are associated with this cell. */
    short int nr_mm_tasks;

  } grav;

  /*! Stars variables */
  struct {

    /*! Pointer to the #spart data. */
    struct spart *parts;

    /*! Nr of #spart in this cell. */
    int count;

    /*! Dependency implicit task for the star ghost  (in->ghost->out)*/
    struct task *ghost_in;

    /*! Dependency implicit task for the star ghost  (in->ghost->out)*/
    struct task *ghost_out;

    /*! The star ghost task itself */
    struct task *ghost;

    /*! Linked list of the tasks computing this cell's star density. */
    struct link *density;

    /*! Number of #spart updated in this cell. */
    int updated;

    /*! Is the #spart data of this cell being used in a sub-cell? */
    int hold;

    /*! Spin lock for various uses (#spart case). */
    swift_lock_type lock;

  } stars;

#ifdef WITH_MPI
  /*! MPI variables */
  struct {

    struct {
      /* Task receiving hydro data (positions). */
      struct task *recv_xv;

      /* Task receiving hydro data (density). */
      struct task *recv_rho;

      /* Task receiving hydro data (gradient). */
      struct task *recv_gradient;

      /* Linked list for sending hydro data (positions). */
      struct link *send_xv;

      /* Linked list for sending hydro data (density). */
      struct link *send_rho;

      /* Linked list for sending hydro data (gradient). */
      struct link *send_gradient;

    } hydro;

    struct {

      /* Task receiving gpart data. */
      struct task *recv;

      /* Linked list for sending gpart data. */
      struct link *send;
    } grav;

    /* Task receiving data (time-step). */
    struct task *recv_ti;

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

  } mpi;
#endif

  /*! The first kick task */
  struct task *kick1;

  /*! The second kick task */
  struct task *kick2;

  /*! The task to compute time-steps */
  struct task *timestep;

  /*! Task for source terms */
  struct task *sourceterms;

  /*! Minimum dimension, i.e. smallest edge of this cell (min(width)). */
  float dmin;

  /*! ID of the previous owner, e.g. runner. */
  int owner;

  /*! ID of the node this cell lives on. */
  int nodeID;

  /*! Number of tasks that are associated with this cell. */
  short int nr_tasks;

  /*! The depth of this cell in the tree. */
  char depth;

  /*! Is this cell split ? */
  char split;

  /*! The maximal depth of this cell and its progenies */
  char maxdepth;

#ifdef SWIFT_DEBUG_CHECKS
  /* Cell ID (for debugging) */
  int cellID;

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
int cell_pack(struct cell *c, struct pcell *pc, const int with_gravity);
int cell_unpack(struct pcell *pc, struct cell *c, struct space *s,
                const int with_gravity);
int cell_pack_tags(const struct cell *c, int *tags);
int cell_unpack_tags(const int *tags, struct cell *c);
int cell_pack_end_step(struct cell *c, struct pcell_step *pcell);
int cell_unpack_end_step(struct cell *c, struct pcell_step *pcell);
int cell_pack_multipoles(struct cell *c, struct gravity_tensors *m);
int cell_unpack_multipoles(struct cell *c, struct gravity_tensors *m);
int cell_getsize(struct cell *c);
int cell_link_parts(struct cell *c, struct part *parts);
int cell_link_gparts(struct cell *c, struct gpart *gparts);
int cell_link_sparts(struct cell *c, struct spart *sparts);
void cell_clean_links(struct cell *c, void *data);
void cell_make_multipoles(struct cell *c, integertime_t ti_current);
void cell_check_multipole(struct cell *c);
void cell_check_foreign_multipole(const struct cell *c);
void cell_clean(struct cell *c);
void cell_check_part_drift_point(struct cell *c, void *data);
void cell_check_gpart_drift_point(struct cell *c, void *data);
void cell_check_multipole_drift_point(struct cell *c, void *data);
void cell_reset_task_counters(struct cell *c);
int cell_unskip_hydro_tasks(struct cell *c, struct scheduler *s);
int cell_unskip_stars_tasks(struct cell *c, struct scheduler *s);
int cell_unskip_gravity_tasks(struct cell *c, struct scheduler *s);
void cell_set_super(struct cell *c, struct cell *super);
void cell_drift_part(struct cell *c, const struct engine *e, int force);
void cell_drift_gpart(struct cell *c, const struct engine *e, int force);
void cell_drift_multipole(struct cell *c, const struct engine *e);
void cell_drift_all_multipoles(struct cell *c, const struct engine *e);
void cell_check_timesteps(struct cell *c);
void cell_store_pre_drift_values(struct cell *c);
void cell_activate_subcell_hydro_tasks(struct cell *ci, struct cell *cj,
                                       struct scheduler *s);
void cell_activate_subcell_grav_tasks(struct cell *ci, struct cell *cj,
                                      struct scheduler *s);
void cell_activate_subcell_stars_tasks(struct cell *ci, struct cell *cj,
                                       struct scheduler *s);
void cell_activate_subcell_external_grav_tasks(struct cell *ci,
                                               struct scheduler *s);
void cell_activate_drift_part(struct cell *c, struct scheduler *s);
void cell_activate_drift_gpart(struct cell *c, struct scheduler *s);
void cell_activate_sorts(struct cell *c, int sid, struct scheduler *s);
void cell_clear_drift_flags(struct cell *c, void *data);
void cell_set_super_mapper(void *map_data, int num_elements, void *extra_data);
int cell_has_tasks(struct cell *c);
void cell_remove_part(const struct engine *e, struct cell *c, struct part *p,
                      struct xpart *xp);
void cell_remove_gpart(const struct engine *e, struct cell *c,
                       struct gpart *gp);
void cell_remove_spart(const struct engine *e, struct cell *c,
                       struct spart *sp);
void cell_convert_part_to_gpart(const struct engine *e, struct cell *c,
                                struct part *p, struct xpart *xp);
void cell_convert_spart_to_gpart(const struct engine *e, struct cell *c,
                                 struct spart *sp);
int cell_can_use_pair_mm(const struct cell *ci, const struct cell *cj,
                         const struct engine *e, const struct space *s);
int cell_can_use_pair_mm_rebuild(const struct cell *ci, const struct cell *cj,
                                 const struct engine *e, const struct space *s);

/* Inlined functions (for speed). */

/**
 * @brief Can a sub-pair hydro task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_recurse_in_pair_hydro_task(const struct cell *c) {

  /* Is the cell split ? */
  /* If so, is the cut-off radius plus the max distance the parts have moved */
  /* smaller than the sub-cell sizes ? */
  /* Note: We use the _old values as these might have been updated by a drift */
  return c->split && ((kernel_gamma * c->hydro.h_max_old +
                       c->hydro.dx_max_part_old) < 0.5f * c->dmin);
}

/**
 * @brief Can a sub-self hydro task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_recurse_in_self_hydro_task(const struct cell *c) {

  /* Is the cell split and not smaller than the smoothing length? */
  return c->split && (kernel_gamma * c->hydro.h_max_old < 0.5f * c->dmin);
}

/**
 * @brief Can a sub-pair star task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_recurse_in_pair_stars_task(const struct cell *c) {

  // LOIC: To implement
  return 0;
}

/**
 * @brief Can a sub-self stars task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_recurse_in_self_stars_task(const struct cell *c) {

  // LOIC: To implement
  return 0;
}

/**
 * @brief Can a pair hydro task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_can_split_pair_hydro_task(
    const struct cell *c) {

  /* Is the cell split ? */
  /* If so, is the cut-off radius with some leeway smaller than */
  /* the sub-cell sizes ? */
  /* Note that since tasks are create after a rebuild no need to take */
  /* into account any part motion (i.e. dx_max == 0 here) */
  return c->split &&
         (space_stretch * kernel_gamma * c->hydro.h_max < 0.5f * c->dmin);
}

/**
 * @brief Can a self hydro task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_can_split_self_hydro_task(
    const struct cell *c) {

  /* Is the cell split ? */
  /* If so, is the cut-off radius with some leeway smaller than */
  /* the sub-cell sizes ? */
  /* Note: No need for more checks here as all the sub-pairs and sub-self */
  /* tasks will be created. So no need to check for h_max */
  return c->split &&
         (space_stretch * kernel_gamma * c->hydro.h_max < 0.5f * c->dmin);
}

/**
 * @brief Can a pair stars task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_can_split_pair_stars_task(
    const struct cell *c) {

  // LOIC: To implement
  return 0;
}

/**
 * @brief Can a self stars task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_can_split_self_stars_task(
    const struct cell *c) {

  // LOIC: To implement
  return 0;
}

/**
 * @brief Can a pair gravity task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_split_pair_gravity_task(const struct cell *c) {

  /* Is the cell split ? */
  return c->split && c->depth < space_subdepth_grav;
}

/**
 * @brief Can a self gravity task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_split_self_gravity_task(const struct cell *c) {

  /* Is the cell split ? */
  return c->split && c->depth < space_subdepth_grav;
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
  return (kernel_gamma * max(ci->hydro.h_max, cj->hydro.h_max) +
              ci->hydro.dx_max_part + cj->hydro.dx_max_part >
          cj->dmin);
}

/**
 * @brief Add a unique tag to a cell mostly for MPI communications.
 *
 * @param c The #cell to tag.
 */
__attribute__((always_inline)) INLINE static void cell_tag(struct cell *c) {
#ifdef WITH_MPI

#ifdef SWIFT_DEBUG_CHECKS
  if (c->mpi.tag > 0) error("setting tag for already tagged cell");
#endif

  if (c->mpi.tag < 0 &&
      (c->mpi.tag = atomic_inc(&cell_next_tag)) > cell_max_tag)
    error("Ran out of cell tags.");
#else
  error("SWIFT was not compiled with MPI enabled.");
#endif  // WITH_MPI
}

#endif /* SWIFT_CELL_H */
