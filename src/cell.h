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
#ifndef SWIFT_CELL_H
#define SWIFT_CELL_H

/* Config parameters. */
#include <config.h>

/* Includes. */
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* Local includes. */
#include "align.h"
#include "cell_black_holes.h"
#include "cell_grav.h"
#include "cell_hydro.h"
#include "cell_rt.h"
#include "cell_sinks.h"
#include "cell_stars.h"
#include "ghost_stats.h"
#include "kernel_hydro.h"
#include "multipole_struct.h"
#include "part.h"
#include "periodic.h"
#include "sort_part.h"
#include "space.h"
#include "task.h"
#include "timeline.h"
#include "zoom_region/zoom.h"

/* Avoid cyclic inclusions */
struct engine;
struct scheduler;
struct replication_list;

/* Max tag size set to 2^29 to take into account some MPI implementations
 * that use 2^31 as the upper bound on MPI tags and the fact that
 * cell_next_tag is multiplied by 2 when passed to an MPI function.
 * The maximum was lowered by a further factor of 2 to be on the safe side.*/
#define cell_max_tag (1 << 29)

#define cell_align 128

/* Global variables. */
extern int cell_next_tag;

/*! Counter for cell IDs (when exceeding max values for uniqueness) */
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
extern unsigned long long last_cell_id;
extern unsigned long long last_leaf_cell_id;
#endif

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

/* Holds the pairs of progeny for each sid. */
struct cell_split_pair {
  int count;
  struct {
    int pid;
    int pjd;
    int sid;
  } pairs[16];
};
extern struct cell_split_pair cell_split_pairs[13];

/**
 * @brief Packed cell for information correct at rebuild time.
 *
 * Contains all the information for a tree walk in a non-local cell.
 */
struct pcell {

  /*! Hydro variables */
  struct {

    /*! Number of #part in this cell. */
    int count;

    /*! Maximal smoothing length. */
    float h_max;

    /*! Minimal integer end-of-timestep in this cell for hydro tasks */
    integertime_t ti_end_min;

    /*! Maximal integer beginning-of-timestep in this cell for hydro tasks */
    integertime_t ti_beg_max;

    /*! Integer time of the last drift of the #part in this cell */
    integertime_t ti_old_part;

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

    /*! Maximal smoothing length. */
    float h_max;

    /*! Minimal integer end-of-timestep in this cell for stars tasks */
    integertime_t ti_end_min;

    /*! Integer time of the last drift of the #spart in this cell */
    integertime_t ti_old_part;

  } stars;

  /*! Black hole variables */
  struct {

    /*! Number of #spart in this cell. */
    int count;

    /*! Maximal smoothing length. */
    float h_max;

    /*! Minimal integer end-of-timestep in this cell for black hole tasks */
    integertime_t ti_end_min;

    /*! Integer time of the last drift of the #spart in this cell */
    integertime_t ti_old_part;

  } black_holes;

  /*! Sink variables */
  struct {

    /*! Number of #sink in this cell. */
    int count;

    /*! Maximal cut off radius. */
    float r_cut_max;

    /*! Minimal integer end-of-timestep in this cell for sinks tasks */
    integertime_t ti_end_min;

    /*! Integer time of the last drift of the #sink in this cell */
    integertime_t ti_old_part;

  } sinks;

  /*! RT variables */
  struct {

    /*! Minimal integer end-of-timestep in this cell for RT tasks */
    integertime_t ti_rt_end_min;

    /*! smallest RT time-step size in this cell */
    integertime_t ti_rt_min_step_size;

  } rt;

  /*! Maximal depth in that part of the tree */
  int maxdepth;

  /*! Relative indices of the cell's progeny. */
  int progeny[8];

#ifdef SWIFT_DEBUG_CHECKS
  /* Cell ID (for debugging) */
  unsigned long long cellID;
#endif

} SWIFT_STRUCT_ALIGN;

/**
 * @brief Cell information at the end of a time-step.
 */
struct pcell_step {

  struct {

    /*! Minimal integer end-of-timestep in this cell (hydro) */
    integertime_t ti_end_min;

    /*! Maximal distance any #part has travelled since last rebuild */
    float dx_max_part;
  } hydro;

  struct {

    /*! Minimal integer end-of-timestep in this cell (gravity) */
    integertime_t ti_end_min;
  } grav;

  struct {

    /*! Minimal integer end-of-timestep in this cell (stars) */
    integertime_t ti_end_min;

    /*! Maximal distance any #part has travelled since last rebuild */
    float dx_max_part;
  } stars;

  struct {

    /*! Minimal integer end-of-timestep in this cell (black_holes) */
    integertime_t ti_end_min;

    /*! Maximal distance any #part has travelled since last rebuild */
    float dx_max_part;
  } black_holes;

  struct {

    /*! Minimal integer end-of-timestep in this cell (rt) */
    integertime_t ti_rt_end_min;

    /*! smallest RT time-step size in this cell */
    integertime_t ti_rt_min_step_size;

  } rt;
};

/**
 * @brief Cell information to propagate the new counts of star particles.
 */
struct pcell_sf {

  /*! Stars variables */
  struct {

    /* Distance by which the stars pointer has moved since the last rebuild */
    ptrdiff_t delta_from_rebuild;

    /* Number of particles in the cell */
    int count;

    /*! Maximum part movement in this cell since last construction. */
    float dx_max_part;

  } stars;

  /*! Grav. variables */
  struct {

    /* Distance by which the gpart pointer has moved since the last rebuild */
    ptrdiff_t delta_from_rebuild;

    /* Number of particles in the cell */
    int count;

  } grav;
};

/**
 * @brief Bitmasks for the cell flags. Beware when adding flags that you don't
 * exceed the size of the flags variable in the struct cell.
 */
enum cell_flags {
  cell_flag_split = (1UL << 0),
  cell_flag_do_hydro_drift = (1UL << 1),
  cell_flag_do_hydro_sub_drift = (1UL << 2),
  cell_flag_do_hydro_sub_sort = (1UL << 3),
  cell_flag_do_hydro_limiter = (1UL << 4),
  cell_flag_do_hydro_sub_limiter = (1UL << 5),
  cell_flag_do_grav_drift = (1UL << 6),
  cell_flag_do_grav_sub_drift = (1UL << 7),
  cell_flag_do_stars_sub_sort = (1UL << 8),
  cell_flag_do_stars_drift = (1UL << 9),
  cell_flag_do_stars_sub_drift = (1UL << 10),
  cell_flag_do_bh_drift = (1UL << 11),
  cell_flag_do_bh_sub_drift = (1UL << 12),
  cell_flag_do_sink_drift = (1UL << 13),
  cell_flag_do_sink_sub_drift = (1UL << 14),
  cell_flag_do_stars_resort = (1UL << 15),
  cell_flag_has_tasks = (1UL << 16),
  cell_flag_do_hydro_sync = (1UL << 17),
  cell_flag_do_hydro_sub_sync = (1UL << 18),
  cell_flag_unskip_self_grav_processed = (1UL << 19),
  cell_flag_unskip_pair_grav_processed = (1UL << 20),
  cell_flag_skip_rt_sort = (1UL << 21),    /* skip rt_sort after a RT recv? */
  cell_flag_do_rt_sub_sort = (1UL << 22),  /* same as hydro_sub_sort for RT */
  cell_flag_rt_requests_sort = (1UL << 23) /* was this sort requested by RT? */
};

/**
 * @brief Names of the cell types.
 */
extern const char *cellID_names[];

/**
 * @brief Names of the cell sub-types.
 */
extern const char *subcellID_names[];

/**
 * @brief What type of top level cell is this cell?
 *
 * All cells are cell_type_regular when running a periodic box. The other types
 * are never used in a periodic box and conversely, cell_type_regular cells are
 * never used when running with a zoom region.
 *
 * When running with a zoom region:
 *
 * - Background cells are low resolution cells covering the majority of
 *   the volume. The zoom region fills a number of these cells in the centre
 *   of the volume (when buffer cells are not used, see below).
 * - Buffer cells are only used when explicitly turned on by the user. These
 *   are a high resolution type of background cell (but are lower resolution
 *   than zoom cells) that are used to pad the volume between the zoom region
 *   and the background cells containing the zoom region.
 * - Zoom cells are the high resolution cells that cover the zoom region
 *   (nested inside the central background/buffer cell/s).
 */
enum cell_types {
  cell_type_regular, /* A standard top level cell (for non-zoom boxes). */
  cell_type_zoom,    /* A zoom cell (only applicable for zooms). */
  cell_type_buffer,  /* A buffer cell (only applicable for zooms). */
  cell_type_bkg,     /* A background cell (only applicable for zooms). */
} __attribute__((__packed__));

/**
 * @brief What subtype of top level cell is this cell?
 *
 * When running a periodic box, cells can only have regular_sub type.
 *
 * When running with a zoom region:
 *
 * - Zoom cells can only be cell_subtype_regular.
 * - Buffer cells (if turned on) can be neighbours if they are within the
 *   gravity criterion of the zoom region or void cells if they contain the
 *   zoom region. Otherwise, they are cell_subtype_regular.
 * - Like buffer cells, background cells can be neighbours if they are within
 *   the gravity criterion of the zoom region or void cells if they contain the
 *   zoom region, but only if buffer cells are not turned on. If buffer cells
 *   are turned on, a background cell can be cell_subtype_empty if it contains
 *   nested buffer cells (only background cells can be empty). Otherwise, they
 *   are cell_subtype_regular.
 *
 * All cell types serve a function but only cell_subtype_neighbour and
 * cell_subtype_regular can get tasks.
 *
 * Void cells do not contain any pointers to particles but carry multipoles and
 * particle counts based on the nested zoom cells.
 *
 * Empty cells do not contain anything, they should not feature in any
 * calculation and only exist to ensure the cell grids are maintained.
 */
enum cell_subtypes {
  cell_subtype_regular,   /* A normal cell. */
  cell_subtype_neighbour, /* A cell within the gravity criterion of the zoom
                             region. */
  cell_subtype_void,      /* A cell containing the zoom region (void cell). */
  cell_subtype_empty      /* An empty cell (background cells containing buffer
                             cells). */
} __attribute__((__packed__));

/**
 * @brief Cell within the tree structure.
 *
 * Contains particles, links to tasks, a multipole object and counters.
 */
struct cell {

  /*! The cell location on the grid (corner nearest to the origin). */
  double loc[3];

  /*! The cell dimensions. */
  double width[3];

  /*! Pointers to the next level of cells. */
  struct cell *progeny[8];

  union {

    /*! Linking pointer for "memory management". */
    struct cell *next;

    /*! Parent cell. */
    struct cell *parent;
  };

  /*! Pointer to the top-level cell in a hierarchy */
  struct cell *top;

  /*! Super cell, i.e. the highest-level parent cell with *any* task */
  struct cell *super;

  /*! The direct void cell parent of a zoom cell. Only used if running with
   * a zoom region. */
  struct cell *void_parent;

  /*! Cell flags bit-mask. */
  volatile uint32_t flags;

  /*! Hydro variables */
  struct cell_hydro hydro;

  /*! Grav variables */
  struct cell_grav grav;

  /*! Stars variables */
  struct cell_stars stars;

  /*! Black hole variables */
  struct cell_black_holes black_holes;

  /*! Sink particles variables */
  struct cell_sinks sinks;

  /*! Radiative transfer variables */
  struct cell_rt rt;

#ifdef WITH_MPI
  /*! MPI variables */
  struct {

    union {
      /* Single list of all send tasks associated with this cell. */
      struct link *send;

      /* Single list of all recv tasks associated with this cell. */
      struct link *recv;
    };

    union {
      /* Single list of all pack tasks associated with this cell. */
      struct link *pack;

      /* Single list of all unpack tasks associated with this cell. */
      struct link *unpack;
    };

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

  /*! The task to limit the time-step of inactive particles */
  struct task *timestep_limiter;

  /*! The task to synchronize the time-step of inactive particles hit by
   * feedback */
  struct task *timestep_sync;

  /*! The task to recursively collect time-steps */
  struct task *timestep_collect;

#ifdef WITH_CSDS
  /*! The csds task */
  struct task *csds;
#endif

  /*! Minimum dimension, i.e. smallest edge of this cell (min(width)). */
  float dmin;

  /*! ID of the previous owner, e.g. runner. */
  short int owner;

  /*! ID of a threadpool thread that maybe associated with this cell. */
  short int tpid;

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

  /*! What type of cell is this? */
  enum cell_types type;

  /*! What subtype of cell is this? */
  enum cell_subtypes subtype;

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
  /* Cell ID (for debugging) */
  long long cellID;
#endif

#ifdef SWIFT_DEBUG_CHECKS

  /*! The list of tasks that have been executed on this cell */
  char tasks_executed[task_type_count];

  /*! The list of sub-tasks that have been executed on this cell */
  char subtasks_executed[task_type_count];
#endif

  struct ghost_stats ghost_statistics;

} SWIFT_STRUCT_ALIGN;

/* Convert cell location to ID. */
#define cell_getid(cdim, i, j, k) \
  ((int)(k) + (cdim)[2] * ((int)(j) + (cdim)[1] * (int)(i)))

/* Function prototypes. */
void cell_split(struct cell *c, ptrdiff_t parts_offset, ptrdiff_t sparts_offset,
                ptrdiff_t bparts_offset, ptrdiff_t sinks_offset,
                struct cell_buff *buff, struct cell_buff *sbuff,
                struct cell_buff *bbuff, struct cell_buff *gbuff,
                struct cell_buff *sinkbuff);
void cell_sanitize(struct cell *c, int treated);
int cell_locktree(struct cell *c);
void cell_unlocktree(struct cell *c);
int cell_glocktree(struct cell *c);
void cell_gunlocktree(struct cell *c);
int cell_mlocktree(struct cell *c);
void cell_munlocktree(struct cell *c);
int cell_slocktree(struct cell *c);
void cell_sunlocktree(struct cell *c);
int cell_sink_locktree(struct cell *c);
void cell_sink_unlocktree(struct cell *c);
int cell_blocktree(struct cell *c);
void cell_bunlocktree(struct cell *c);
int cell_pack(struct cell *c, struct pcell *pc, const int with_gravity);
int cell_unpack(struct pcell *pc, struct cell *c, struct space *s,
                const int with_gravity);
void cell_pack_part_swallow(const struct cell *c,
                            struct black_holes_part_data *data);
void cell_unpack_part_swallow(struct cell *c,
                              const struct black_holes_part_data *data);
void cell_pack_bpart_swallow(const struct cell *c,
                             struct black_holes_bpart_data *data);
void cell_unpack_bpart_swallow(struct cell *c,
                               const struct black_holes_bpart_data *data);
int cell_pack_tags(const struct cell *c, int *tags);
int cell_unpack_tags(const int *tags, struct cell *c);
int cell_pack_end_step(const struct cell *c, struct pcell_step *pcell);
int cell_unpack_end_step(struct cell *c, const struct pcell_step *pcell);
void cell_pack_timebin(const struct cell *const c, timebin_t *const t);
void cell_unpack_timebin(struct cell *const c, timebin_t *const t);
int cell_pack_multipoles(struct cell *c, struct gravity_tensors *m);
int cell_unpack_multipoles(struct cell *c, struct gravity_tensors *m);
int cell_pack_sf_counts(struct cell *c, struct pcell_sf *pcell);
int cell_unpack_sf_counts(struct cell *c, struct pcell_sf *pcell);
int cell_get_tree_size(struct cell *c);
int cell_link_parts(struct cell *c, struct part *parts);
int cell_link_gparts(struct cell *c, struct gpart *gparts);
int cell_link_sparts(struct cell *c, struct spart *sparts);
int cell_link_bparts(struct cell *c, struct bpart *bparts);
int cell_link_foreign_parts(struct cell *c, struct part *parts);
int cell_link_foreign_gparts(struct cell *c, struct gpart *gparts);
void cell_unlink_foreign_particles(struct cell *c);
int cell_count_parts_for_tasks(const struct cell *c);
int cell_count_gparts_for_tasks(const struct cell *c);
void cell_clean_links(struct cell *c, void *data);
void cell_make_multipoles(struct cell *c, integertime_t ti_current,
                          const struct gravity_props *const grav_props);
void cell_check_multipole(struct cell *c,
                          const struct gravity_props *const grav_props);
void cell_check_foreign_multipole(const struct cell *c);
void cell_clean(struct cell *c);
void cell_check_part_drift_point(struct cell *c, void *data);
void cell_check_gpart_drift_point(struct cell *c, void *data);
void cell_check_spart_drift_point(struct cell *c, void *data);
void cell_check_sink_drift_point(struct cell *c, void *data);
void cell_check_multipole_drift_point(struct cell *c, void *data);
void cell_reset_task_counters(struct cell *c);
int cell_unskip_hydro_tasks(struct cell *c, struct scheduler *s);
int cell_unskip_stars_tasks(struct cell *c, struct scheduler *s,
                            const int with_star_formation,
                            const int with_star_formation_sink);
int cell_unskip_sinks_tasks(struct cell *c, struct scheduler *s);
int cell_unskip_rt_tasks(struct cell *c, struct scheduler *s,
                         const int sub_cycle);
int cell_unskip_black_holes_tasks(struct cell *c, struct scheduler *s);
int cell_unskip_gravity_tasks(struct cell *c, struct scheduler *s);
void cell_drift_part(struct cell *c, const struct engine *e, int force,
                     struct replication_list *replication_list_in);
void cell_drift_gpart(struct cell *c, const struct engine *e, int force,
                      struct replication_list *replication_list);
void cell_drift_spart(struct cell *c, const struct engine *e, int force,
                      struct replication_list *replication_list);
void cell_drift_sink(struct cell *c, const struct engine *e, int force);
void cell_drift_bpart(struct cell *c, const struct engine *e, int force,
                      struct replication_list *replication_list);
void cell_drift_multipole(struct cell *c, const struct engine *e);
void cell_drift_all_multipoles(struct cell *c, const struct engine *e);
void cell_check_timesteps(const struct cell *c, const integertime_t ti_current,
                          const timebin_t max_bin);
void cell_store_pre_drift_values(struct cell *c);
void cell_set_star_resort_flag(struct cell *c);
void cell_activate_star_formation_tasks(struct cell *c, struct scheduler *s,
                                        const int with_feedback);
void cell_activate_star_formation_sink_tasks(struct cell *c,
                                             struct scheduler *s,
                                             const int with_feedback);
void cell_activate_sink_formation_tasks(struct cell *c, struct scheduler *s);
void cell_activate_subcell_hydro_tasks(struct cell *ci, struct cell *cj,
                                       struct scheduler *s,
                                       const int with_timestep_limiter);
int cell_activate_subcell_grav_tasks(struct cell *ci, struct cell *cj,
                                     struct scheduler *s);
void cell_activate_subcell_stars_tasks(struct cell *ci, struct cell *cj,
                                       struct scheduler *s,
                                       const int with_star_formation,
                                       const int with_star_formation_sink,
                                       const int with_timestep_sync);
void cell_activate_subcell_sinks_tasks(struct cell *ci, struct cell *cj,
                                       struct scheduler *s,
                                       const int with_timestep_sync);
void cell_activate_subcell_black_holes_tasks(struct cell *ci, struct cell *cj,
                                             struct scheduler *s,
                                             const int with_timestep_sync);
void cell_activate_subcell_external_grav_tasks(struct cell *ci,
                                               struct scheduler *s);
void cell_activate_subcell_rt_tasks(struct cell *ci, struct cell *cj,
                                    struct scheduler *s, const int sub_cycle);
void cell_set_no_rt_sort_flag_up(struct cell *c);
void cell_activate_super_spart_drifts(struct cell *c, struct scheduler *s);
void cell_activate_super_sink_drifts(struct cell *c, struct scheduler *s);
void cell_activate_drift_part(struct cell *c, struct scheduler *s);
void cell_activate_drift_gpart(struct cell *c, struct scheduler *s);
void cell_activate_drift_spart(struct cell *c, struct scheduler *s);
void cell_activate_drift_sink(struct cell *c, struct scheduler *s);
void cell_activate_drift_bpart(struct cell *c, struct scheduler *s);
void cell_activate_sync_part(struct cell *c, struct scheduler *s);
void cell_activate_rt_sorts(struct cell *c, int sid, struct scheduler *s);
void cell_activate_hydro_sorts(struct cell *c, int sid, struct scheduler *s);
void cell_activate_stars_sorts(struct cell *c, int sid, struct scheduler *s);
void cell_activate_limiter(struct cell *c, struct scheduler *s);
void cell_clear_drift_flags(struct cell *c, void *data);
void cell_clear_limiter_flags(struct cell *c, void *data);
void cell_set_super_mapper(void *map_data, int num_elements, void *extra_data);
void cell_check_spart_pos(const struct cell *c,
                          const struct spart *global_sparts);
void cell_check_sort_flags(const struct cell *c);
void cell_clear_stars_sort_flags(struct cell *c, const int unused_flags);
void cell_clear_hydro_sort_flags(struct cell *c, const int unused_flags);
void cell_clear_unskip_flags(struct cell *c);
int cell_has_tasks(struct cell *c);
void cell_remove_part(const struct engine *e, struct cell *c, struct part *p,
                      struct xpart *xp);
void cell_remove_gpart(const struct engine *e, struct cell *c,
                       struct gpart *gp);
void cell_remove_spart(const struct engine *e, struct cell *c,
                       struct spart *sp);
void cell_remove_sink(const struct engine *e, struct cell *c,
                      struct sink *sink);
void cell_remove_bpart(const struct engine *e, struct cell *c,
                       struct bpart *bp);
struct spart *cell_add_spart(struct engine *e, struct cell *c);
struct gpart *cell_add_gpart(struct engine *e, struct cell *c);
struct spart *cell_spawn_new_spart_from_part(struct engine *e, struct cell *c,
                                             const struct part *p,
                                             const struct xpart *xp);
struct spart *cell_spawn_new_spart_from_sink(struct engine *e, struct cell *c,
                                             const struct sink *s);
struct gpart *cell_convert_part_to_gpart(const struct engine *e, struct cell *c,
                                         struct part *p, struct xpart *xp);
struct gpart *cell_convert_spart_to_gpart(const struct engine *e,
                                          struct cell *c, struct spart *sp);
struct gpart *cell_convert_bpart_to_gpart(const struct engine *e,
                                          struct cell *c, struct bpart *bp);
struct spart *cell_convert_part_to_spart(struct engine *e, struct cell *c,
                                         struct part *p, struct xpart *xp);
struct sink *cell_convert_part_to_sink(struct engine *e, struct cell *c,
                                       struct part *p, struct xpart *xp);
void cell_reorder_extra_parts(struct cell *c, const ptrdiff_t parts_offset);
void cell_reorder_extra_gparts(struct cell *c, struct part *parts,
                               struct spart *sparts, struct sink *sinks);
void cell_reorder_extra_sparts(struct cell *c, const ptrdiff_t sparts_offset);
void cell_reorder_extra_sinks(struct cell *c, const ptrdiff_t sinks_offset);
int cell_can_use_pair_mm(const struct cell *ci, const struct cell *cj,
                         const struct engine *e, const struct space *s,
                         const int use_rebuild_data, const int is_tree_walk,
                         const int periodic, const int use_mesh);

/***
 * @brief Get the cell ID of a cell including an offset.
 *
 * NOTE: This function is only used in the zoom code.
 *
 * @param cdim The dimensions of the cell grid.
 * @param offset The offset to add to the cell ID.
 * @param i, j, k The cell ijk coordinates.
 *
 * @return The cell id.
 * */
__attribute__((always_inline)) INLINE int cell_getid_offset(const int cdim[3],
                                                            const int offset,
                                                            const int i,
                                                            const int j,
                                                            const int k) {
  return cell_getid(cdim, i, j, k) + offset;
}

/**
 * @brief For a given location, what TL cell does it belong to in one of the
 * cell grids below the background level?
 *
 * NOTE: This function is only applicable when running with a zoom region.
 *
 * It finds the cell ID in a region below the background level (zoom/buffer)
 * using the lower boundary of the nested region to offset the input coordinates
 * before multiplying by the inverse width of a cell to get the integer cell
 * coordinates.
 *
 * This function is mainly just a convenience wrapper around cell_getid_offset
 * to remove some boilerplate when applying the bounds.
 *
 * @param cdim The cell grid dimensions.
 * @param bounds The edges of this nested region.
 * @param x, y, z Location of particle within buffer or zoom region.
 * @param iwidth The width of a cell in this grid.
 * @param offset The offset of this cell type in cells_top.
 *
 * @return The cell id.
 */
__attribute__((always_inline)) INLINE int cell_getid_below_bkg(
    const int cdim[3], const double bounds[3], const double x, const double y,
    const double z, const double iwidth[3], const int offset) {

  /* Get the cell ijk coordinates in this grid. */
  const int i = (x - bounds[0]) * iwidth[0];
  const int j = (y - bounds[1]) * iwidth[1];
  const int k = (z - bounds[2]) * iwidth[2];

  /* Which zoom TL cell are we in? */
  return cell_getid_offset(cdim, offset, i, j, k);
}

/**
 * @brief For a given particle location, what TL cell does it belong to?
 *
 * Any calls to cell_getid_from_pos will be redirected to this function when
 * running with a zoom region. This function will then identify which level of
 * top level (TL) cell should be returned, either background, buffer (if used),
 * or zoom.
 *
 * The cell hierarchy is structured with background cells at the top (largest),
 * followed by buffer cells (intermediate), and finally zoom cells (smallest).
 * When running without buffer cells the hierarchy is bkg -> zoom.
 *
 * We do this by testing each level from the top down. First we see what
 * background cell the position is in. If it's a void cell, we then get the
 * zoom cell (in the no buffer cell case). If it's an empty cell (buffer cell
 * case), we then get the buffer cell. If it's then a void buffer cell, we
 * then get the zoom region.
 *
 * Cells are not stored in their hierarchy order in s->cells_top. Instead,
 * zoom cells are first, followed by background cells, and finally buffer cells.
 * This is a bit strange but puts the zoom cells as the primary cell type and
 * means we can simplify the code surrounding hydro operations which are
 * isolated to the zoom region.
 *
 * @param s The space.
 * @param x, y, z Location to get the cell ID for.
 *
 * @return The cell id.
 */
__attribute__((always_inline)) INLINE int zoom_cell_getid(const struct space *s,
                                                          const double x,
                                                          const double y,
                                                          const double z) {

  /* Lets get some properties of the zoom region. */
  const struct zoom_region_properties *zoom_props = s->zoom_props;
  const int bkg_cell_offset = zoom_props->bkg_cell_offset;
  const double zoom_lower_bounds[3] = {zoom_props->region_lower_bounds[0],
                                       zoom_props->region_lower_bounds[1],
                                       zoom_props->region_lower_bounds[2]};
  const int buffer_cell_offset = zoom_props->buffer_cell_offset;
  const double buffer_lower_bounds[3] = {zoom_props->buffer_lower_bounds[0],
                                         zoom_props->buffer_lower_bounds[1],
                                         zoom_props->buffer_lower_bounds[2]};

  /* Get the background cell ijk coordinates. */
  const int bkg_i = x * s->iwidth[0];
  const int bkg_j = y * s->iwidth[1];
  const int bkg_k = z * s->iwidth[2];

  /* Which background cell is this? */
  int cell_id =
      cell_getid_offset(s->cdim, bkg_cell_offset, bkg_i, bkg_j, bkg_k);

#ifdef SWIFT_DEBUG_CHECKS
  if (cell_id < 0 || cell_id >= s->nr_cells)
    error("cell_id out of range: %i (%f %f %f)", cell_id, x, y, z);
#endif

  /* If this is a void cell we are in the zoom region. */
  if (s->cells_top[cell_id].subtype == cell_subtype_void) {

    /* Which zoom TL cell are we in? */
    return cell_getid_below_bkg(s->zoom_props->cdim, zoom_lower_bounds, x, y, z,
                                s->zoom_props->iwidth,
                                /*offset*/ 0);

  }

  /* If this is an empty cell we are in the buffer cells.
   * Otherwise, It's a legitimate background cell, and we'll return it. */
  else if (s->cells_top[cell_id].subtype == cell_subtype_empty) {

    /* Which buffer TL cell are we in? */
    cell_id = cell_getid_below_bkg(
        s->zoom_props->buffer_cdim, buffer_lower_bounds, x, y, z,
        s->zoom_props->buffer_iwidth, buffer_cell_offset);

    /* Here we need to check if this is the void buffer cell.
     * Otherwise, It's a legitimate buffer cell, and we'll return it. */
    if (s->cells_top[cell_id].subtype == cell_subtype_void) {

      /* Which zoom TL cell are we in? */
      return cell_getid_below_bkg(s->zoom_props->cdim, zoom_lower_bounds, x, y,
                                  z, s->zoom_props->iwidth,
                                  /*offset*/ 0);
    }
  }

  return cell_id;
}

/**
 * @brief Convert a coordinate to the cell ID containing it.
 *
 * @param s The space.
 * @param x, y, z Coordinates of particle/cell.
 *
 * @return The cell id.
 */
__attribute__((always_inline)) INLINE int cell_getid_from_pos(
    const struct space *s, const double x, const double y, const double z) {

  /* When running with a zoom region we need to account for the nested
   * cell grids, so call the zoom specific version. */
  if (s->with_zoom_region) {

    /* Use the version that accounts for the zoom region */
    return zoom_cell_getid(s, x, y, z);

  } else {

    /* Zoom region isn't enabled so we can use the simple version */
    const int i = x * s->iwidth[0];
    const int j = y * s->iwidth[1];
    const int k = z * s->iwidth[2];
    return cell_getid(s->cdim, i, j, k);
  }
}

/**
 * @brief Does a #cell contain no particle at all.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_is_empty(
    const struct cell *c) {

  /* Void cells are never empty. */
  if (c->subtype == cell_subtype_void) return 0;

  return (c->hydro.count == 0 && c->grav.count == 0 && c->stars.count == 0 &&
          c->black_holes.count == 0 && c->sinks.count == 0);
}

/**
 * @brief Compute the square of the distance between the CoMs of two multipoles.
 *
 * @param multi_i The first #gravity_tensors.
 * @param multi_j The second #gravity_tensors.
 * @param use_rebuild_data Are we considering the data at the last tree-build
 * (1) or current data (0)?
 * @param periodic Are we using periodic BCs?
 * @param dim The dimensions of the simulation volume.
 *
 * @return The square of the distance between the multiepoles' CoMs.
 */
__attribute__((always_inline)) INLINE static double cell_mpole_CoM_dist2(
    const struct gravity_tensors *restrict multi_i,
    const struct gravity_tensors *restrict multi_j, const int use_rebuild_data,
    const int periodic, const double dim[3]) {

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
  return dx * dx + dy * dy + dz * dz;
}

/**
 * @brief Compute the square of the minimal distance between any two points in
 * two cells of the same size
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param periodic Are we using periodic BCs?
 * @param dim The dimensions of the simulation volume
 *
 * @return The square of the minimal distance between the two cells.
 */
__attribute__((always_inline)) INLINE static double cell_min_dist2(
    const struct cell *restrict ci, const struct cell *restrict cj,
    const int periodic, const double dim[3]) {

  const double cix_min = ci->loc[0];
  const double ciy_min = ci->loc[1];
  const double ciz_min = ci->loc[2];
  const double cjx_min = cj->loc[0];
  const double cjy_min = cj->loc[1];
  const double cjz_min = cj->loc[2];

  const double cix_max = ci->loc[0] + ci->width[0];
  const double ciy_max = ci->loc[1] + ci->width[1];
  const double ciz_max = ci->loc[2] + ci->width[2];
  const double cjx_max = cj->loc[0] + cj->width[0];
  const double cjy_max = cj->loc[1] + cj->width[1];
  const double cjz_max = cj->loc[2] + cj->width[2];

  if (periodic) {

    const double dx = min4(fabs(nearest(cix_min - cjx_min, dim[0])),
                           fabs(nearest(cix_min - cjx_max, dim[0])),
                           fabs(nearest(cix_max - cjx_min, dim[0])),
                           fabs(nearest(cix_max - cjx_max, dim[0])));

    const double dy = min4(fabs(nearest(ciy_min - cjy_min, dim[1])),
                           fabs(nearest(ciy_min - cjy_max, dim[1])),
                           fabs(nearest(ciy_max - cjy_min, dim[1])),
                           fabs(nearest(ciy_max - cjy_max, dim[1])));

    const double dz = min4(fabs(nearest(ciz_min - cjz_min, dim[2])),
                           fabs(nearest(ciz_min - cjz_max, dim[2])),
                           fabs(nearest(ciz_max - cjz_min, dim[2])),
                           fabs(nearest(ciz_max - cjz_max, dim[2])));

    return dx * dx + dy * dy + dz * dz;

  } else {

    const double dx = min4(fabs(cix_min - cjx_min), fabs(cix_min - cjx_max),
                           fabs(cix_max - cjx_min), fabs(cix_max - cjx_max));
    const double dy = min4(fabs(ciy_min - cjy_min), fabs(ciy_min - cjy_max),
                           fabs(ciy_max - cjy_min), fabs(ciy_max - cjy_max));
    const double dz = min4(fabs(ciz_min - cjz_min), fabs(ciz_min - cjz_max),
                           fabs(ciz_max - cjz_min), fabs(ciz_max - cjz_max));

    return dx * dx + dy * dy + dz * dz;
  }
}

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
 * @param ci The #cell with stars.
 * @param cj The #cell with hydro parts.
 */
__attribute__((always_inline)) INLINE static int
cell_can_recurse_in_pair_stars_task(const struct cell *ci,
                                    const struct cell *cj) {

  /* Is the cell split ? */
  /* If so, is the cut-off radius plus the max distance the parts have moved */
  /* smaller than the sub-cell sizes ? */
  /* Note: We use the _old values as these might have been updated by a drift */
  return ci->split && cj->split &&
         ((kernel_gamma * ci->stars.h_max_old + ci->stars.dx_max_part_old) <
          0.5f * ci->dmin) &&
         ((kernel_gamma * cj->hydro.h_max_old + cj->hydro.dx_max_part_old) <
          0.5f * cj->dmin);
}

/**
 * @brief Can a sub-self stars task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_recurse_in_self_stars_task(const struct cell *c) {

  /* Is the cell split and not smaller than the smoothing length? */
  return c->split && (kernel_gamma * c->stars.h_max_old < 0.5f * c->dmin) &&
         (kernel_gamma * c->hydro.h_max_old < 0.5f * c->dmin);
}

/**
 * @brief Can a sub-pair sink task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param ci The #cell with stars.
 * @param cj The #cell with hydro parts.
 */
__attribute__((always_inline)) INLINE static int
cell_can_recurse_in_pair_sinks_task(const struct cell *ci,
                                    const struct cell *cj) {

  /* Is the cell split ? */
  /* If so, is the cut-off radius plus the max distance the parts have moved */
  /* smaller than the sub-cell sizes ? */
  /* Note: We use the _old values as these might have been updated by a drift */
  return ci->split && cj->split &&
         ((ci->sinks.r_cut_max_old + ci->sinks.dx_max_part_old) <
          0.5f * ci->dmin) &&
         ((kernel_gamma * cj->hydro.h_max_old + cj->hydro.dx_max_part_old) <
          0.5f * cj->dmin);
}

/**
 * @brief Can a sub-pair black hole task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param ci The #cell with black holes.
 * @param cj The #cell with hydro parts.
 */
__attribute__((always_inline)) INLINE static int
cell_can_recurse_in_pair_black_holes_task(const struct cell *ci,
                                          const struct cell *cj) {

  /* Is the cell split ? */
  /* If so, is the cut-off radius plus the max distance the parts have moved */
  /* smaller than the sub-cell sizes ? */
  /* Note: We use the _old values as these might have been updated by a drift */
  return ci->split && cj->split &&
         ((kernel_gamma * ci->black_holes.h_max_old +
           ci->black_holes.dx_max_part_old) < 0.5f * ci->dmin) &&
         ((kernel_gamma * cj->hydro.h_max_old + cj->hydro.dx_max_part_old) <
          0.5f * cj->dmin);
}

/**
 * @brief Can a sub-self black hole task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_recurse_in_self_black_holes_task(const struct cell *c) {

  /* Is the cell split and not smaller than the smoothing length? */
  return c->split &&
         (kernel_gamma * c->black_holes.h_max_old < 0.5f * c->dmin) &&
         (kernel_gamma * c->hydro.h_max_old < 0.5f * c->dmin);
}

/**
 * @brief Can a sub-self sinks task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_recurse_in_self_sinks_task(const struct cell *c) {

  /* Is the cell split and not smaller than the smoothing length? */
  return c->split && (c->sinks.r_cut_max_old < 0.5f * c->dmin) &&
         (kernel_gamma * c->hydro.h_max_old < 0.5f * c->dmin);
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
         (space_stretch * kernel_gamma * c->hydro.h_max < 0.5f * c->dmin) &&
         (space_stretch * kernel_gamma * c->stars.h_max < 0.5f * c->dmin) &&
         (space_stretch * kernel_gamma * c->black_holes.h_max < 0.5f * c->dmin);
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
         (space_stretch * kernel_gamma * c->hydro.h_max < 0.5f * c->dmin) &&
         (space_stretch * kernel_gamma * c->stars.h_max < 0.5f * c->dmin) &&
         (space_stretch * kernel_gamma * c->black_holes.h_max < 0.5f * c->dmin);
}

/**
 * @brief Is a cell above the task level?
 *
 * The task level is defined as space_subdepth_diff_grav above the leaves
 * of the tree.
 *
 * When running a zoom simulation this will use zoom_bkg_subdepth_diff_grav for
 * the background cells while the zoom cells will use the regular threshold.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_is_above_diff_grav_depth(
    const struct cell *c) {

  /* Void cells must always be treated all the ways to the leaves so we will
   * always return true here. */
  if (c->subtype == cell_subtype_void) {
    return 1;
  }

  /* Regular and zoom cells use the usual condition. */
  if (c->type == cell_type_regular || c->type == cell_type_zoom) {
    return (c->maxdepth - c->depth) > space_subdepth_diff_grav;
  }

  /* Otherwise all other cells use the background diff_grav constant. */
  return (c->maxdepth - c->depth) > zoom_bkg_subdepth_diff_grav;
}

/**
 * @brief Can a pair gravity task associated with a pair of cells be split
 * into smaller sub-tasks?
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_split_pair_gravity_task(const struct cell *ci, const struct cell *cj) {

  /* Otherwise, are the cells split and still far from the leaves ? */
  return (ci->split && cj->split) && cell_is_above_diff_grav_depth(ci) &&
         cell_is_above_diff_grav_depth(cj);
}

/**
 * @brief Can a pair gravity task be split into smaller sub-tasks
 * based on the number of particles in the cells?
 *
 * If the product of the number of particles (the number of interactions) is
 * smaller than space_subsize_pair_grav a task will not be split since it is
 * already small enough and this function will return 1.
 *
 * Void cells will always be split regardless of the number of particles.
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_pair_gravity_task_below_subsize(const struct cell *ci,
                                     const struct cell *cj) {

  /* Get the cell counts. */
  const long long gcount_i = ci->grav.count;
  const long long gcount_j = cj->grav.count;

  return gcount_i * gcount_j < ((long long)space_subsize_pair_grav);
}

/**
 * @brief Can a self gravity task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_split_self_gravity_task(const struct cell *c) {

  /* Is the cell split and still far from the leaves ? */
  return c->split && cell_is_above_diff_grav_depth(c);
}

/**
 * @brief Can a self FOF task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_can_split_self_fof_task(
    const struct cell *c) {

  /* Is the cell split ? */
  return c->split && c->grav.count > 5000 && cell_is_above_diff_grav_depth(c);
}

/**
 * @brief Have gas particles in a pair of cells moved too much and require a
 * rebuild
 * ?
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
__attribute__((always_inline, nonnull)) INLINE static int
cell_need_rebuild_for_hydro_pair(const struct cell *ci, const struct cell *cj) {

  /* Is the cut-off radius plus the max distance the parts in both cells have */
  /* moved larger than the cell size ? */
  /* Note ci->dmin == cj->dmin */
  if (kernel_gamma * max(ci->hydro.h_max, cj->hydro.h_max) +
          ci->hydro.dx_max_part + cj->hydro.dx_max_part >
      cj->dmin) {
    return 1;
  }
  return 0;
}

/**
 * @brief Have star particles in a pair of cells moved too much and require a
 * rebuild?
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
__attribute__((always_inline, nonnull)) INLINE static int
cell_need_rebuild_for_stars_pair(const struct cell *ci, const struct cell *cj) {

  /* Is the cut-off radius plus the max distance the parts in both cells have */
  /* moved larger than the cell size ? */
  /* Note ci->dmin == cj->dmin */
  if (kernel_gamma * max(ci->stars.h_max, cj->hydro.h_max) +
          ci->stars.dx_max_part + cj->hydro.dx_max_part >
      cj->dmin) {
    return 1;
  }
  return 0;
}

/**
 * @brief Have sink particles in a pair of cells moved too much and require a
 * rebuild?
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
__attribute__((always_inline, nonnull)) INLINE static int
cell_need_rebuild_for_sinks_pair(const struct cell *ci, const struct cell *cj) {

  /* Is the cut-off radius plus the max distance the parts in both cells have */
  /* moved larger than the cell size ? */
  /* Note ci->dmin == cj->dmin */
  if (max(ci->sinks.r_cut_max, kernel_gamma * cj->hydro.h_max) +
          ci->sinks.dx_max_part + cj->hydro.dx_max_part >
      cj->dmin) {
    return 1;
  }
  return 0;
}

/**
 * @brief Have star particles in a pair of cells moved too much and require a
 * rebuild?
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
__attribute__((always_inline, nonnull)) INLINE static int
cell_need_rebuild_for_black_holes_pair(const struct cell *ci,
                                       const struct cell *cj) {

  /* Is the cut-off radius plus the max distance the parts in both cells have */
  /* moved larger than the cell size ? */
  /* Note ci->dmin == cj->dmin */
  if (kernel_gamma * max(ci->black_holes.h_max, cj->hydro.h_max) +
          ci->black_holes.dx_max_part + cj->hydro.dx_max_part >
      cj->dmin) {
    return 1;
  }
  return 0;
}

/**
 * @brief Add a unique tag to a cell, mostly for MPI communications.
 *
 * This function locks the cell so that tags can be added concurrently.
 *
 * @param c The #cell to tag.
 */
__attribute__((always_inline)) INLINE static void cell_ensure_tagged(
    struct cell *c) {
#ifdef WITH_MPI

  lock_lock(&c->hydro.lock);
  if (c->mpi.tag < 0 &&
      (c->mpi.tag = atomic_inc(&cell_next_tag)) > cell_max_tag)
    error("Ran out of cell tags.");
  if (lock_unlock(&c->hydro.lock) != 0) {
    error("Failed to unlock cell.");
  }
#else
  error("SWIFT was not compiled with MPI enabled.");
#endif  // WITH_MPI
}

/**
 * @brief Allocate hydro sort memory for cell.
 *
 * @param c The #cell that will require sorting.
 * @param flags Cell flags.
 */
__attribute__((always_inline)) INLINE static void cell_malloc_hydro_sorts(
    struct cell *c, const int flags) {

  const int count = c->hydro.count;

  /* Have we already allocated something? */
  if (c->hydro.sort != NULL) {

    /* Start by counting how many dimensions we need
       and how many we already have */
    const int num_arrays_wanted =
        intrinsics_popcount(c->hydro.sort_allocated | flags);
    const int num_already_allocated =
        intrinsics_popcount(c->hydro.sort_allocated);

    /* Do we already have what we want? */
    if (num_arrays_wanted == num_already_allocated) return;

    /* Allocate memory for the new array */
    struct sort_entry *new_array = NULL;
    if ((new_array = (struct sort_entry *)swift_malloc(
             "hydro.sort", sizeof(struct sort_entry) * num_arrays_wanted *
                               (count + 1))) == NULL)
      error("Failed to allocate sort memory.");

    /* Now, copy the already existing arrays */
    int from = 0;
    int to = 0;
    for (int j = 0; j < 13; j++) {
      if (c->hydro.sort_allocated & (1 << j)) {
        memcpy(&new_array[to * (count + 1)], &c->hydro.sort[from * (count + 1)],
               sizeof(struct sort_entry) * (count + 1));
        ++from;
        ++to;
      } else if (flags & (1 << j)) {
        ++to;
        c->hydro.sort_allocated |= (1 << j);
      }
    }

    /* Swap the pointers */
    swift_free("hydro.sort", c->hydro.sort);
    c->hydro.sort = new_array;

  } else {

    c->hydro.sort_allocated = flags;

    /* Start by counting how many dimensions we need */
    const int num_arrays = intrinsics_popcount(flags);

    /* If there is anything, allocate enough memory */
    if (num_arrays) {
      if ((c->hydro.sort = (struct sort_entry *)swift_malloc(
               "hydro.sort",
               sizeof(struct sort_entry) * num_arrays * (count + 1))) == NULL)
        error("Failed to allocate sort memory.");
    }
  }
}

/**
 * @brief Free hydro sort memory for cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static void cell_free_hydro_sorts(
    struct cell *c) {

#ifdef NONE_SPH
  /* Nothing to do as we have no particles and hence no sorts */
#else
  if (c->hydro.sort != NULL) {
    swift_free("hydro.sort", c->hydro.sort);
    c->hydro.sort = NULL;
    c->hydro.sort_allocated = 0;
  }
#endif
}

/**
 * @brief Returns the array of sorted indices for the gas particles of a given
 * cell along agiven direction.
 *
 * @param c The #cell.
 * @param sid the direction id.
 */
__attribute__((always_inline)) INLINE static struct sort_entry *
cell_get_hydro_sorts(const struct cell *c, const int sid) {

#ifdef SWIFT_DEBUG_CHECKS
  if (sid >= 13 || sid < 0) error("Invalid sid!");

  if (!(c->hydro.sort_allocated & (1 << sid)))
    error("Sort not allocated along direction %d", sid);
#endif

  /* We need to find at what position in the meta-array of
     sorts where the corresponding sid has been allocated since
     there might be gaps as we only allocated the directions that
     are in use.
     We create a mask with all the bits before the sid's one set to 1
     and apply it on the list of allocated directions. We then count
     the number of bits that are in the results to obtain the position
     of the correspondin sid in the meta-array */
  const int j = intrinsics_popcount(c->hydro.sort_allocated & ((1 << sid) - 1));

  /* Return the corresponding array */
  return &c->hydro.sort[j * (c->hydro.count + 1)];
}

/**
 * @brief Allocate stars sort memory for cell.
 *
 * @param c The #cell that will require sorting.
 * @param flags Cell flags.
 */
__attribute__((always_inline)) INLINE static void cell_malloc_stars_sorts(
    struct cell *c, const int flags) {

  const int count = c->stars.count;

  /* Have we already allocated something? */
  if (c->stars.sort != NULL) {

    /* Start by counting how many dimensions we need
       and how many we already have */
    const int num_arrays_wanted =
        intrinsics_popcount(c->stars.sort_allocated | flags);
    const int num_already_allocated =
        intrinsics_popcount(c->stars.sort_allocated);

    /* Do we already have what we want? */
    if (num_arrays_wanted == num_already_allocated) return;

    /* Allocate memory for the new array */
    struct sort_entry *new_array = NULL;
    if ((new_array = (struct sort_entry *)swift_malloc(
             "stars.sort", sizeof(struct sort_entry) * num_arrays_wanted *
                               (count + 1))) == NULL)
      error("Failed to allocate sort memory.");

    /* Now, copy the already existing arrays */
    int from = 0;
    int to = 0;
    for (int j = 0; j < 13; j++) {
      if (c->stars.sort_allocated & (1 << j)) {
        memcpy(&new_array[to * (count + 1)], &c->stars.sort[from * (count + 1)],
               sizeof(struct sort_entry) * (count + 1));
        ++from;
        ++to;
      } else if (flags & (1 << j)) {
        ++to;
        c->stars.sort_allocated |= (1 << j);
      }
    }

    /* Swap the pointers */
    swift_free("stars.sort", c->stars.sort);
    c->stars.sort = new_array;

  } else {

    c->stars.sort_allocated = flags;

    /* Start by counting how many dimensions we need */
    const int num_arrays = intrinsics_popcount(flags);

    /* If there is anything, allocate enough memory */
    if (num_arrays) {
      if ((c->stars.sort = (struct sort_entry *)swift_malloc(
               "stars.sort",
               sizeof(struct sort_entry) * num_arrays * (count + 1))) == NULL)
        error("Failed to allocate sort memory.");
    }
  }
}

/**
 * @brief Free stars sort memory for cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static void cell_free_stars_sorts(
    struct cell *c) {

#ifdef STARS_NONE
  /* Nothing to do as we have no particles and hence no sorts */
#else
  if (c->stars.sort != NULL) {
    swift_free("stars.sort", c->stars.sort);
    c->stars.sort = NULL;
    c->stars.sort_allocated = 0;
  }
#endif
}

/**
 * @brief Returns the array of sorted indices for the star particles of a given
 * cell along agiven direction.
 *
 * @param c The #cell.
 * @param sid the direction id.
 */
__attribute__((always_inline)) INLINE static struct sort_entry *
cell_get_stars_sorts(const struct cell *c, const int sid) {

#ifdef SWIFT_DEBUG_CHECKS
  if (sid >= 13 || sid < 0) error("Invalid sid!");

  if (!(c->stars.sort_allocated & (1 << sid)))
    error("Sort not allocated along direction %d", sid);
#endif

  /* We need to find at what position in the meta-array of
     sorts where the corresponding sid has been allocated since
     there might be gaps as we only allocated the directions that
     are in use.
     We create a mask with all the bits before the sid's one set to 1
     and apply it on the list of allocated directions. We then count
     the number of bits that are in the results to obtain the position
     of the correspondin sid in the meta-array */
  const int j = intrinsics_popcount(c->stars.sort_allocated & ((1 << sid) - 1));

  /* Return the corresponding array */
  return &c->stars.sort[j * (c->stars.count + 1)];
}

/**
 * @brief Set the given flag for the given cell.
 */
__attribute__((always_inline)) INLINE static void cell_set_flag(
    struct cell *c, const uint32_t flag) {
  atomic_or(&c->flags, flag);
}

/**
 * @brief Clear the given flag for the given cell.
 */
__attribute__((always_inline)) INLINE static void cell_clear_flag(
    struct cell *c, const uint32_t flag) {
  atomic_and(&c->flags, ~flag);
}

/**
 * @brief  Get the given flag for the given cell.
 */
__attribute__((always_inline)) INLINE static int cell_get_flag(
    const struct cell *c, const uint32_t flag) {
  return (c->flags & flag) > 0;
}

/**
 * @brief Check if a cell has a recv task of the given subtype.
 */
__attribute__((always_inline)) INLINE static struct task *cell_get_recv(
    const struct cell *c, enum task_subtypes subtype) {
#ifdef WITH_MPI
  struct link *l = c->mpi.recv;
  while (l != NULL && l->t->subtype != subtype) l = l->next;
  return (l != NULL) ? l->t : NULL;
#else
  return NULL;
#endif
}

/**
 * @brief Generate the cell ID for top level cells. Only used for debugging.
 *
 * Cell IDs are stored in the long long `cell->cellID`. Top level cells get
 * their index according to their location on the top level grid.
 * We have 15 bits set aside in `cell->cellID` for the top level cells. Hence
 * if we have more that 32^3 top level cells, the cell IDs won't be guaranteed
 * to be unique and reproducible between two runs, but only unique.
 *
 * @param c #cell to work with
 * @param s The #space
 */
__attribute__((always_inline)) INLINE void cell_assign_top_level_cell_index(
    struct cell *c, struct space *s) {

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
  if (c->depth != 0) {
    error("assigning top level cell index to cell with depth > 0");
  } else {

    if (!s->with_zoom_region) {

      /* Unpack properties we'll need */
      const int *cdim = s->cdim;
      const double *iwidth = s->iwidth;

      if (cdim[0] * cdim[1] * cdim[2] > 32 * 32 * 32) {
        /* print warning only once */
        if (last_cell_id == 1ULL) {
          message(
              "WARNING: Got %d x %d x %d top level cells. "
              "Cell IDs are only guaranteed to be "
              "reproduceably unique if count is < 32^3",
              cdim[0], cdim[1], cdim[2]);
        }
        /* Do this in same line. Otherwise, bad things happen. */
        c->cellID = atomic_inc(&last_cell_id);
      } else {
        int i = (int)(c->loc[0] * iwidth[0] + 0.5);
        int j = (int)(c->loc[1] * iwidth[1] + 0.5);
        int k = (int)(c->loc[2] * iwidth[2] + 0.5);
        c->cellID = (unsigned long long)(cell_getid(cdim, i, j, k) + 1);
      }
      /* in both cases, keep track of first prodigy index */
      atomic_inc(&last_leaf_cell_id);
    } else {

      /* In the zoom case we assign 2 or 3 separate cell grids (3 when buffer
       * cells have been turned on, 2 otherwise). These are all stored in
       * s->cells_top with zoom cells first, followed by background cells, and
       * finally buffer cells (if used).
       *
       * Therefore: zoom_cell_ids < background_cell_ids < buffer_cell_ids
       *
       * The offsets of each cell grid are stored for accessing the cells, but
       * this means cell ids are only reproducable if the total number of all
       * cells is less than 32^3. */

      /* Get the cdims */
      const int cdim[3] = {s->cdim[0], s->cdim[1], s->cdim[2]};
      const int zoom_cdim[3] = {s->zoom_props->cdim[0], s->zoom_props->cdim[1],
                                s->zoom_props->cdim[2]};
      const int buffer_cdim[3] = {s->zoom_props->buffer_cdim[0],
                                  s->zoom_props->buffer_cdim[1],
                                  s->zoom_props->buffer_cdim[2]};

      if (((cdim[0] * cdim[1] * cdim[2]) +
           (zoom_cdim[0] * zoom_cdim[1] * zoom_cdim[2]) +
           (buffer_cdim[0] * buffer_cdim[1] * buffer_cdim[2])) > 32 * 32 * 32) {
        /* print warning only once */
        if (last_cell_id == 1ULL) {
          message(
              "WARNING: Got (%d x %d x %d + %d x %d x %d + %d x %d x %d) top "
              "level cells. "
              "Cell IDs are only guaranteed to be "
              "reproduceably unique if count is < 32^3",
              cdim[0], cdim[1], cdim[2], zoom_cdim[0], zoom_cdim[1],
              zoom_cdim[2], buffer_cdim[0], buffer_cdim[1], buffer_cdim[2]);
        }
        /* Do this in same line. Otherwise, bad things happen. */
        c->cellID = atomic_inc(&last_cell_id);
      } else {
        c->cellID = (unsigned long long)(cell_getid_from_pos(
            s, c->loc[0], c->loc[1], c->loc[2]));
      }
      /* in both cases, keep track of first prodigy index */
      atomic_inc(&last_leaf_cell_id);
    }
  }
#endif
}

/**
 * @brief Generate the cell ID for progeny cells. Only used for debugging.
 *
 * Cell IDs are stored in the unsigned long long `cell->cellID`.
 * We have 15 bits set aside in `cell->cellID` for the top level cells, with
 * 49 remaining. Each progeny cell gets a unique ID by inheriting
 * its parent ID and adding 3 bits on the left side, which are set according
 * to the progeny's location within its parent cell. Finally, a 1 is set as the
 * leading bit such that all recursive children with index (000) are still
 * recognized as such. This allows us to give IDs to 16 levels of depth
 * uniquely.
 * If the depth exceeds 16, we use the old scheme where we just add up a
 * counter, which is not a reproducible way of giving IDs to cells, but
 * guarantees uniqueness.
 */
__attribute__((always_inline)) INLINE void cell_assign_cell_index(
    struct cell *c, const struct cell *parent) {

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
  if (c->depth == 0) error("assigning progeny cell index to top level cell.");
  if (c->depth > 16 || last_cell_id > 1ULL) {
    /* last_cell_id > 1 => too many top level cells for clever IDs */
    if (last_cell_id == 1ULL) { /* warning not yet printed; do it only once */
      message(
          "WARNING: Got depth %d > 16."
          "IDs are only guaranteed unique if depth <= 16",
          c->depth);
      last_cell_id += 1ULL;
    }
    /* Do this in same line. Otherwise, bad things happen. */
    c->cellID = atomic_inc(&last_leaf_cell_id);
  } else {
    /* we're good to go for unique IDs */
    /* first inherit the parent's ID */
    unsigned long long child_id = parent->cellID;

    /* if parent isn't top level cell, we have to
     * remove the marker (leading 1) of the previous depth first,
     * as we're going to add 3 bits for this new depth at that
     * position in the variable now. So turn that leading 1 into a 0 */
    if (c->depth > 1) child_id &= ~(1ULL << ((c->depth - 1) * 3 + 15));

    /* Now add marker (leading 1) for this depth 3 bits further to the left*/
    child_id |= 1ULL << (15 + c->depth * 3);

    /* get progeny index in parent cell */
    unsigned long long local_id = 0LL;
    if (c->loc[0] > parent->loc[0]) local_id |= 1LL;
    if (c->loc[1] > parent->loc[1]) local_id |= 2LL;
    if (c->loc[2] > parent->loc[2]) local_id |= 4LL;
    local_id <<= (15 + (c->depth - 1) * 3);

    /* add progeny index to cell index */
    child_id |= local_id;

    c->cellID = child_id;
  }

#endif
}

#endif /* SWIFT_CELL_H */
