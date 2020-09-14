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
#include <stdint.h>
#include <string.h>

/* Local includes. */
#include "align.h"
#include "kernel_hydro.h"
#include "kernel_dark_matter.h"
#include "lock.h"
#include "multipole_struct.h"
#include "part.h"
#include "periodic.h"
#include "sort_part.h"
#include "space.h"
#include "star_formation_logger_struct.h"
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

    /*! Maximal integer end-of-timestep in this cell for hydro tasks */
    integertime_t ti_end_max;

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

    /*! Maximal smoothing length. */
    float h_max;

    /*! Minimal integer end-of-timestep in this cell for stars tasks */
    integertime_t ti_end_min;

    /*! Maximal integer end-of-timestep in this cell for stars tasks */
    integertime_t ti_end_max;

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

    /*! Maximal integer end-of-timestep in this cell for black hole tasks */
    integertime_t ti_end_max;

    /*! Integer time of the last drift of the #spart in this cell */
    integertime_t ti_old_part;

  } black_holes;
    
    /*! Dark matter variables */
    struct {
        
        /*! Number of #spart in this cell. */
        int count;
        
        /*! Maximal smoothing length. */
        float h_max;
        
        /*! Minimal integer end-of-timestep in this cell for stars tasks */
        integertime_t ti_end_min;
        
        /*! Maximal integer end-of-timestep in this cell for stars tasks */
        integertime_t ti_end_max;
        
        /*! Integer time of the last drift of the #spart in this cell */
        integertime_t ti_old_part;
        
    } dark_matter;

  /*! Maximal depth in that part of the tree */
  int maxdepth;

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
struct pcell_step_hydro {

  /*! Minimal integer end-of-timestep in this cell (hydro) */
  integertime_t ti_end_min;

  /*! Minimal integer end-of-timestep in this cell (hydro) */
  integertime_t ti_end_max;

  /*! Maximal distance any #part has travelled since last rebuild */
  float dx_max_part;
};

struct pcell_step_grav {

  /*! Minimal integer end-of-timestep in this cell (gravity) */
  integertime_t ti_end_min;

  /*! Minimal integer end-of-timestep in this cell (gravity) */
  integertime_t ti_end_max;
};

struct pcell_step_stars {

  /*! Minimal integer end-of-timestep in this cell (stars) */
  integertime_t ti_end_min;

  /*! Maximal integer end-of-timestep in this cell (stars) */
  integertime_t ti_end_max;

  /*! Maximal distance any #part has travelled since last rebuild */
  float dx_max_part;
};

struct pcell_step_black_holes {

  /*! Minimal integer end-of-timestep in this cell (black_holes) */
  integertime_t ti_end_min;

  /*! Maximal integer end-of-timestep in this cell (black_holes) */
  integertime_t ti_end_max;

  /*! Maximal distance any #part has travelled since last rebuild */
  float dx_max_part;
};

struct pcell_step_dark_matter {
    
    /*! Minimal integer end-of-timestep in this cell (dark_matter) */
    integertime_t ti_end_min;
    
    /*! Maximal integer end-of-timestep in this cell (dark_matter) */
    integertime_t ti_end_max;
    
    /*! Maximal distance any #part has travelled since last rebuild */
    float dx_max_part;
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
  cell_flag_do_stars_resort = (1UL << 13),
  cell_flag_has_tasks = (1UL << 14),
  cell_flag_do_hydro_sync = (1UL << 15),
  cell_flag_do_hydro_sub_sync = (1UL << 16),
  cell_flag_do_dark_matter_drift = (1UL << 17),
  cell_flag_do_dark_matter_sub_drift = (1UL << 18),
  cell_flag_do_dark_matter_sub_sort = (1UL << 19)
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

  /*! Pointers to the next level of cells. */
  struct cell *progeny[8];

  /*! Linking pointer for "memory management". */
  struct cell *next;

  /*! Parent cell. */
  struct cell *parent;

  /*! Pointer to the top-level cell in a hierarchy */
  struct cell *top;

  /*! Super cell, i.e. the highest-level parent cell with *any* task */
  struct cell *super;

  /*! Cell flags bit-mask. */
  volatile uint32_t flags;

  /*! Hydro variables */
  struct {

    /*! Pointer to the #part data. */
    struct part *parts;

    /*! Pointer to the #xpart data. */
    struct xpart *xparts;

    /*! Pointer for the sorted indices. */
    struct sort_entry *sort;

    /*! Super cell, i.e. the highest-level parent cell that has a hydro
     * pair/self tasks */
    struct cell *super;

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

    /*! Linked list of the tasks computing this cell's limiter. */
    struct link *limiter;

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

    /*! Dependency implicit task for cooling (in->cooling->out) */
    struct task *cooling_in;

    /*! Dependency implicit task for cooling (in->cooling->out) */
    struct task *cooling_out;

    /*! Task for cooling */
    struct task *cooling;

    /*! Task for star formation */
    struct task *star_formation;

    /*! Task for sorting the stars again after a SF event */
    struct task *stars_resort;

    /*! Last (integer) time the cell's part were drifted forward in time. */
    integertime_t ti_old_part;

    /*! Minimum end of (integer) time step in this cell for hydro tasks. */
    integertime_t ti_end_min;

    /*! Maximum end of (integer) time step in this cell for hydro tasks. */
    integertime_t ti_end_max;

    /*! Maximum beginning of (integer) time step in this cell for hydro tasks.
     */
    integertime_t ti_beg_max;

    /*! Spin lock for various uses (#part case). */
    swift_lock_type lock;

    /*! Max smoothing length in this cell. */
    float h_max;

    /*! Maximum part movement in this cell since last construction. */
    float dx_max_part;

    /*! Maximum particle movement in this cell since the last sort. */
    float dx_max_sort;

    /*! Values of h_max before the drifts, used for sub-cell tasks. */
    float h_max_old;

    /*! Values of dx_max before the drifts, used for sub-cell tasks. */
    float dx_max_part_old;

    /*! Values of dx_max_sort before the drifts, used for sub-cell tasks. */
    float dx_max_sort_old;

    /*! Nr of #part in this cell. */
    int count;

    /*! Nr of #part this cell can hold after addition of new #part. */
    int count_total;

    /*! Number of #part updated in this cell. */
    int updated;

    /*! Is the #part data of this cell being used in a sub-cell? */
    int hold;

    /*! Bit mask of sort directions that will be needed in the next timestep. */
    uint16_t requires_sorts;

    /*! Bit mask of sorts that need to be computed for this cell. */
    uint16_t do_sort;

    /*! Bit-mask indicating the sorted directions */
    uint16_t sorted;

    /*! Bit-mask indicating the sorted directions */
    uint16_t sort_allocated;

#ifdef SWIFT_DEBUG_CHECKS

    /*! Last (integer) time the cell's sort arrays were updated. */
    integertime_t ti_sort;

#endif

  } hydro;

  /*! Grav variables */
  struct {

    /*! Pointer to the #gpart data. */
    struct gpart *parts;

    /*! Pointer to the #spart data at rebuild time. */
    struct gpart *parts_rebuild;

    /*! This cell's multipole. */
    struct gravity_tensors *multipole;

    /*! Super cell, i.e. the highest-level parent cell that has a grav pair/self
     * tasks */
    struct cell *super;
    
    /*! The drift task for gparts */
    struct task *drift;

    /*! Implicit task (going up- and down the tree) for the #gpart drifts */
    struct task *drift_out;

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

    /*! The task to end the force calculation */
    struct task *end_force;

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

    /*! Spin lock for various uses (#gpart case). */
    swift_lock_type plock;

    /*! Spin lock for various uses (#multipole case). */
    swift_lock_type mlock;

    /*! Spin lock for star formation use. */
    swift_lock_type star_formation_lock;

    /*! Nr of #gpart in this cell. */
    int count;

    /*! Nr of #gpart this cell can hold after addition of new #gpart. */
    int count_total;

    /*! Number of #gpart updated in this cell. */
    int updated;

    /*! Is the #gpart data of this cell being used in a sub-cell? */
    int phold;

    /*! Is the #multipole data of this cell being used in a sub-cell? */
    int mhold;

    /*! Number of M-M tasks that are associated with this cell. */
    short int nr_mm_tasks;

  } grav;

  /*! Stars variables */
  struct {

    /*! Pointer to the #spart data. */
    struct spart *parts;

    /*! Pointer to the #spart data at rebuild time. */
    struct spart *parts_rebuild;

    /*! The star ghost task itself */
    struct task *ghost;

    /*! Linked list of the tasks computing this cell's star density. */
    struct link *density;

    /*! Linked list of the tasks computing this cell's star feedback. */
    struct link *feedback;

    /*! The task computing this cell's sorts before the density. */
    struct task *sorts;

    /*! The drift task for sparts */
    struct task *drift;

    /*! Implicit tasks marking the entry of the stellar physics block of tasks
     */
    struct task *stars_in;

    /*! Implicit tasks marking the exit of the stellar physics block of tasks */
    struct task *stars_out;

    /*! Last (integer) time the cell's spart were drifted forward in time. */
    integertime_t ti_old_part;

    /*! Spin lock for various uses (#spart case). */
    swift_lock_type lock;

    /*! Spin lock for star formation use. */
    swift_lock_type star_formation_lock;

    /*! Nr of #spart in this cell. */
    int count;

    /*! Nr of #spart this cell can hold after addition of new #spart. */
    int count_total;

    /*! Max smoothing length in this cell. */
    float h_max;

    /*! Values of h_max before the drifts, used for sub-cell tasks. */
    float h_max_old;

    /*! Maximum part movement in this cell since last construction. */
    float dx_max_part;

    /*! Values of dx_max before the drifts, used for sub-cell tasks. */
    float dx_max_part_old;

    /*! Maximum particle movement in this cell since the last sort. */
    float dx_max_sort;

    /*! Values of dx_max_sort before the drifts, used for sub-cell tasks. */
    float dx_max_sort_old;

    /*! Pointer for the sorted indices. */
    struct sort_entry *sort;

    /*! Bit mask of sort directions that will be needed in the next timestep. */
    uint16_t requires_sorts;

    /*! Bit-mask indicating the sorted directions */
    uint16_t sorted;

    /*! Bit-mask indicating the sorted directions */
    uint16_t sort_allocated;

    /*! Bit mask of sorts that need to be computed for this cell. */
    uint16_t do_sort;

    /*! Maximum end of (integer) time step in this cell for star tasks. */
    integertime_t ti_end_min;

    /*! Maximum end of (integer) time step in this cell for star tasks. */
    integertime_t ti_end_max;

    /*! Maximum beginning of (integer) time step in this cell for star tasks.
     */
    integertime_t ti_beg_max;

    /*! Number of #spart updated in this cell. */
    int updated;

    /*! Is the #spart data of this cell being used in a sub-cell? */
    int hold;

    /*! Star formation history struct */
    struct star_formation_history sfh;

#ifdef SWIFT_DEBUG_CHECKS
    /*! Last (integer) time the cell's sort arrays were updated. */
    integertime_t ti_sort;
#endif

  } stars;
    
    /*! Dark matter variables */
    struct {
        
        /*! Pointer to the #dmpart data. */
        struct dmpart *parts;
        
        /*! Pointer to the #dmpart data at rebuild time. */
        struct dmpart *parts_rebuild;
        
        /*! Pointer for the sorted indices. */
        struct sort_entry *sort;

        /*! The task computing this cell's sorts. */
        struct task *sorts;
        
        /*! Linked list of the tasks computing this cell's dm self-interactions. */
        struct link *sidm;
        
        /*! The drift task for parts */
        struct task *drift;

        /*! The dark matter ghost task itself */
        struct task *ghost;
        
        /*! Linked list of the tasks computing this cell's dark matter density. */
        struct link *density;
        
        /*! kick due to DM-DM interactions */
        struct task *sidm_kick;
        
        /*! Last (integer) time the cell's spart were drifted forward in time. */
        integertime_t ti_old_part;
        
        /*! Spin lock for various uses (#dmpart case). */
        swift_lock_type lock;
        
        /*! Nr of #dmpart in this cell. */
        int count;
        
        /*! Nr of #dmpart this cell can hold after addition of new #dmpart. */
        int count_total;
        
        /*! Max smoothing length in this cell. */
        float h_max;
        
        /*! Values of h_max before the drifts, used for sub-cell tasks. */
        float h_max_old;
        
        /*! Maximum part movement in this cell since last construction. */
        float dx_max_part;
        
        /*! Maximum particle movement in this cell since the last sort. */
        float dx_max_sort;
        
        /*! Values of dx_max before the drifts, used for sub-cell tasks. */
        float dx_max_part_old;
        
        /*! Values of dx_max_sort before the drifts, used for sub-cell tasks. */
        float dx_max_sort_old;
        
        /*! Bit-mask indicating the sorted directions */
        uint16_t sort_allocated;

        /*! Maximum end of (integer) time step in this cell for star tasks. */
        integertime_t ti_end_min;
        
        /*! Maximum end of (integer) time step in this cell for star tasks. */
        integertime_t ti_end_max;
        
        /*! Maximum beginning of (integer) time step in this cell for star tasks.
         */
        integertime_t ti_beg_max;
        
        /*! Number of #spart updated in this cell. */
        int updated;
        
        /*! Is the #dmpart data of this cell being used in a sub-cell? */
        int hold;
        
        /*! Bit mask of sort directions that will be needed in the next timestep. */
        uint16_t requires_sorts;
        
        /*! Bit mask of sorts that need to be computed for this cell. */
        uint16_t do_sort;
        
        /*! Bit-mask indicating the sorted directions */
        uint16_t sorted;
        
        
    } dark_matter;

  /*! Black hole variables */
  struct {

    /*! Pointer to the #bpart data. */
    struct bpart *parts;

    /*! The drift task for bparts */
    struct task *drift;

    /*! Implicit tasks marking the entry of the BH physics block of tasks
     */
    struct task *black_holes_in;

    /*! Implicit tasks marking the exit of the BH physics block of tasks */
    struct task *black_holes_out;

    /*! The star ghost task itself */
    struct task *density_ghost;

    /*! The star ghost task itself */
    struct task *swallow_ghost[3];

    /*! Linked list of the tasks computing this cell's BH density. */
    struct link *density;

    /*! Linked list of the tasks computing this cell's BH swallowing and
     * merging. */
    struct link *swallow;

    /*! Linked list of the tasks processing the particles to swallow */
    struct link *do_gas_swallow;

    /*! Linked list of the tasks processing the particles to swallow */
    struct link *do_bh_swallow;

    /*! Linked list of the tasks computing this cell's BH feedback. */
    struct link *feedback;

    /*! Last (integer) time the cell's bpart were drifted forward in time. */
    integertime_t ti_old_part;

    /*! Spin lock for various uses (#bpart case). */
    swift_lock_type lock;

    /*! Nr of #bpart in this cell. */
    int count;

    /*! Nr of #bpart this cell can hold after addition of new #bpart. */
    int count_total;

    /*! Max smoothing length in this cell. */
    float h_max;

    /*! Values of h_max before the drifts, used for sub-cell tasks. */
    float h_max_old;

    /*! Maximum part movement in this cell since last construction. */
    float dx_max_part;

    /*! Values of dx_max before the drifts, used for sub-cell tasks. */
    float dx_max_part_old;

    /*! Maximum end of (integer) time step in this cell for black tasks. */
    integertime_t ti_end_min;

    /*! Maximum end of (integer) time step in this cell for black hole tasks. */
    integertime_t ti_end_max;

    /*! Maximum beginning of (integer) time step in this cell for black hole
     * tasks.
     */
    integertime_t ti_beg_max;

    /*! Number of #bpart updated in this cell. */
    int updated;

    /*! Is the #bpart data of this cell being used in a sub-cell? */
    int hold;

  } black_holes;

  /*! Sink particles variables */
  struct {

    /*! Pointer to the #sink data. */
    struct sink *parts;

    /*! Nr of #sink in this cell. */
    int count;

    /*! Nr of #sink this cell can hold after addition of new one. */
    int count_total;

    /*! Is the #sink data of this cell being used in a sub-cell? */
    int hold;

    /*! Spin lock for various uses (#sink case). */
    swift_lock_type lock;

    /*! Last (integer) time the cell's sink were drifted forward in time. */
    integertime_t ti_old_part;

    /*! Minimum end of (integer) time step in this cell for sink tasks. */
    integertime_t ti_end_min;

    /*! Maximum end of (integer) time step in this cell for sink tasks. */
    integertime_t ti_end_max;

    /*! Maximum beginning of (integer) time step in this cell for sink
     * tasks.
     */
    integertime_t ti_beg_max;
  } sinks;

#ifdef WITH_MPI
  /*! MPI variables */
  struct {

    union {
      /* Single list of all send tasks associated with this cell. */
      struct link *send;

      /* Single list of all recv tasks associated with this cell. */
      struct link *recv;
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

#ifdef WITH_LOGGER
  /*! The logger task */
  struct task *logger;
#endif

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

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
  /* Cell ID (for debugging) */
  int cellID;
#endif

#ifdef SWIFT_DEBUG_CHECKS

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
                ptrdiff_t bparts_offset, ptrdiff_t dmparts_offset, ptrdiff_t sinks_offset,
                struct cell_buff *buff, struct cell_buff *sbuff,
                struct cell_buff *bbuff, struct cell_buff *gbuff, struct cell_buff *dmbuff,
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
int cell_dmlocktree(struct cell *c);
void cell_dmunlocktree(struct cell *c);
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
int cell_pack_end_step_hydro(struct cell *c, struct pcell_step_hydro *pcell);
int cell_unpack_end_step_hydro(struct cell *c, struct pcell_step_hydro *pcell);
int cell_unpack_end_step_dark_matter(struct cell *c, struct pcell_step_dark_matter *pcell);
int cell_pack_end_step_grav(struct cell *c, struct pcell_step_grav *pcell);
int cell_pack_end_step_dark_matter(struct cell *c, struct pcell_step_dark_matter *pcell);
int cell_unpack_end_step_grav(struct cell *c, struct pcell_step_grav *pcell);
int cell_pack_end_step_stars(struct cell *c, struct pcell_step_stars *pcell);
int cell_unpack_end_step_stars(struct cell *c, struct pcell_step_stars *pcell);
int cell_pack_end_step_black_holes(struct cell *c,
                                   struct pcell_step_black_holes *pcell);
int cell_unpack_end_step_black_holes(struct cell *c,
                                     struct pcell_step_black_holes *pcell);
int cell_pack_multipoles(struct cell *c, struct gravity_tensors *m);
int cell_unpack_multipoles(struct cell *c, struct gravity_tensors *m);
int cell_pack_sf_counts(struct cell *c, struct pcell_sf *pcell);
int cell_unpack_sf_counts(struct cell *c, struct pcell_sf *pcell);
int cell_getsize(struct cell *c);
int cell_link_parts(struct cell *c, struct part *parts);
int cell_link_gparts(struct cell *c, struct gpart *gparts);
int cell_link_sparts(struct cell *c, struct spart *sparts);
int cell_link_bparts(struct cell *c, struct bpart *bparts);
int cell_link_dmparts(struct cell *c, struct dmpart *dmparts);
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
void cell_check_dmpart_drift_point(struct cell *c, void *data);
void cell_check_multipole_drift_point(struct cell *c, void *data);
void cell_reset_task_counters(struct cell *c);
int cell_unskip_hydro_tasks(struct cell *c, struct scheduler *s);
int cell_unskip_stars_tasks(struct cell *c, struct scheduler *s,
                            const int with_star_formation);
int cell_unskip_black_holes_tasks(struct cell *c, struct scheduler *s);
int cell_unskip_dark_matter_tasks(struct cell *c, struct scheduler *s);
int cell_unskip_gravity_tasks(struct cell *c, struct scheduler *s);
void cell_drift_part(struct cell *c, const struct engine *e, int force);
void cell_drift_gpart(struct cell *c, const struct engine *e, int force);
void cell_drift_spart(struct cell *c, const struct engine *e, int force);
void cell_drift_bpart(struct cell *c, const struct engine *e, int force);
void cell_drift_dmpart(struct cell *c, const struct engine *e, int force);
void cell_drift_multipole(struct cell *c, const struct engine *e);
void cell_drift_all_multipoles(struct cell *c, const struct engine *e);
void cell_check_timesteps(const struct cell *c, const integertime_t ti_current,
                          const timebin_t max_bin);
void cell_store_pre_drift_values(struct cell *c);
void cell_set_star_resort_flag(struct cell *c);
void cell_activate_star_formation_tasks(struct cell *c, struct scheduler *s,
                                        const int with_feedback);
void cell_activate_subcell_hydro_tasks(struct cell *ci, struct cell *cj,
                                       struct scheduler *s,
                                       const int with_timestep_limiter);
void cell_activate_subcell_grav_tasks(struct cell *ci, struct cell *cj,
                                      struct scheduler *s);
void cell_activate_subcell_stars_tasks(struct cell *ci, struct cell *cj,
                                       struct scheduler *s,
                                       const int with_star_formation,
                                       const int with_timestep_sync);
void cell_activate_subcell_black_holes_tasks(struct cell *ci, struct cell *cj,
                                             struct scheduler *s,
                                             const int with_timestep_sync);
void cell_activate_subcell_dark_matter_tasks(struct cell *ci, struct cell *cj,
                                             struct scheduler *s);
void cell_activate_subcell_external_grav_tasks(struct cell *ci,
                                               struct scheduler *s);
void cell_activate_super_spart_drifts(struct cell *c, struct scheduler *s);
void cell_activate_drift_part(struct cell *c, struct scheduler *s);
void cell_activate_drift_gpart(struct cell *c, struct scheduler *s);
void cell_activate_drift_spart(struct cell *c, struct scheduler *s);
void cell_activate_drift_bpart(struct cell *c, struct scheduler *s);
void cell_activate_drift_dmpart(struct cell *c, struct scheduler *s);
void cell_activate_sync_part(struct cell *c, struct scheduler *s);
void cell_activate_hydro_sorts(struct cell *c, int sid, struct scheduler *s);
void cell_activate_dark_matter_sorts(struct cell *c, int sid, struct scheduler *s);
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
void cell_clear_dark_matter_sort_flags(struct cell *c, const int unused_flags);
int cell_has_tasks(struct cell *c);
void cell_remove_part(const struct engine *e, struct cell *c, struct part *p,
                      struct xpart *xp);
void cell_remove_gpart(const struct engine *e, struct cell *c,
                       struct gpart *gp);
void cell_remove_spart(const struct engine *e, struct cell *c,
                       struct spart *sp);
void cell_remove_bpart(const struct engine *e, struct cell *c,
                       struct bpart *bp);
struct spart *cell_add_spart(struct engine *e, struct cell *c);
struct gpart *cell_add_gpart(struct engine *e, struct cell *c);
struct spart *cell_spawn_new_spart_from_part(struct engine *e, struct cell *c,
                                             const struct part *p,
                                             const struct xpart *xp);
struct gpart *cell_convert_part_to_gpart(const struct engine *e, struct cell *c,
                                         struct part *p, struct xpart *xp);
struct gpart *cell_convert_spart_to_gpart(const struct engine *e,
                                          struct cell *c, struct spart *sp);
struct spart *cell_convert_part_to_spart(struct engine *e, struct cell *c,
                                         struct part *p, struct xpart *xp);
void cell_reorder_extra_parts(struct cell *c, const ptrdiff_t parts_offset);
void cell_reorder_extra_gparts(struct cell *c, struct part *parts,
                               struct spart *sparts);
void cell_reorder_extra_sparts(struct cell *c, const ptrdiff_t sparts_offset);
int cell_can_use_pair_mm(const struct cell *ci, const struct cell *cj,
                         const struct engine *e, const struct space *s,
                         const int use_rebuild_data, const int is_tree_walk);

/**
 * @brief Compute the square of the minimal distance between any two points in
 * two cells of the same size
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param periodic Are we using periodic BCs?
 * @param dim The dimensions of the simulation volume
 */
__attribute__((always_inline)) INLINE static double cell_min_dist2_same_size(
    const struct cell *restrict ci, const struct cell *restrict cj,
    const int periodic, const double dim[3]) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->width[0] != cj->width[0]) error("Cells of different size!");
  if (ci->width[1] != cj->width[1]) error("Cells of different size!");
  if (ci->width[2] != cj->width[2]) error("Cells of different size!");
#endif

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

    const double dx = min(fabs(cix_max - cjx_min), fabs(cix_min - cjx_max));
    const double dy = min(fabs(ciy_max - cjy_min), fabs(ciy_min - cjy_max));
    const double dz = min(fabs(ciz_max - cjz_min), fabs(ciz_min - cjz_max));

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
 * @brief Can a sub-pair dark_matter task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_recurse_in_pair_dark_matter_task(const struct cell *c) {
    
    /* Is the cell split ? */
    /* If so, is the cut-off radius plus the max distance the parts have moved */
    /* smaller than the sub-cell sizes ? */
    /* Note: We use the _old values as these might have been updated by a drift */
    return c->split && ((dm_kernel_gamma * c->dark_matter.h_max_old +
                         c->dark_matter.dx_max_part_old) < 0.5f * c->dmin);
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
 * @brief Can a sub-self dark matter task recurse to a lower level based
 * on the status of the particles in the cell.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_recurse_in_self_dark_matter_task(const struct cell *c) {
    
    /* Is the cell split and not smaller than the smoothing length? */
    return c->split && (dm_kernel_gamma * c->dark_matter.h_max_old < 0.5f * c->dmin);
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
 * @brief Can a pair hydro task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_can_split_pair_dark_matter_task(const struct cell *c) {
    
    /* Is the cell split ? */
    /* If so, is the cut-off radius with some leeway smaller than */
    /* the sub-cell sizes ? */
    /* Note that since tasks are create after a rebuild no need to take */
    /* into account any part motion (i.e. dx_max == 0 here) */
    return c->split &&
    (space_stretch * dm_kernel_gamma * c->dark_matter.h_max < 0.5f * c->dmin);
}

/**
 * @brief Can a self hydro task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int cell_can_split_self_dark_matter_task(const struct cell *c) {
    
    /* Is the cell split ? */
    /* If so, is the cut-off radius with some leeway smaller than */
    /* the sub-cell sizes ? */
    /* Note: No need for more checks here as all the sub-pairs and sub-self */
    /* tasks will be created. So no need to check for h_max */
    return c->split &&
    (space_stretch * dm_kernel_gamma * c->dark_matter.h_max < 0.5f * c->dmin);
}

/**
 * @brief Can a pair gravity task associated with a cell be split into smaller
 * sub-tasks.
 *
 * @param c The #cell.
 */
__attribute__((always_inline)) INLINE static int
cell_can_split_pair_gravity_task(const struct cell *c) {

  /* Is the cell split and still far from the leaves ? */
  return c->split && ((c->maxdepth - c->depth) > space_subdepth_diff_grav);
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
  return c->split && ((c->maxdepth - c->depth) > space_subdepth_diff_grav);
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
  return c->split && c->grav.count > 5000 &&
         ((c->maxdepth - c->depth) > space_subdepth_diff_grav);
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
 * @brief Have dark matter particles in a pair of cells moved too much and require a
 * rebuild?
 *
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
__attribute__((always_inline, nonnull)) INLINE static int
cell_need_rebuild_for_dark_matter_pair(const struct cell *ci,
                                       const struct cell *cj) {
    
    /* Is the cut-off radius plus the max distance the parts in both cells have */
    /* moved larger than the cell size ? */
    /* Note ci->dmin == cj->dmin */
    if (dm_kernel_gamma * max(ci->dark_matter.h_max, cj->dark_matter.h_max) +
        ci->dark_matter.dx_max_part + cj->dark_matter.dx_max_part >
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

  if (c->hydro.sort != NULL) {
    swift_free("hydro.sort", c->hydro.sort);
    c->hydro.sort = NULL;
    c->hydro.sort_allocated = 0;
  }
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
 * @brief Returns the array of sorted indices for the DM particles of a given
 * cell along agiven direction.
 *
 * @param c The #cell.
 * @param sid the direction id.
 */
__attribute__((always_inline)) INLINE static struct sort_entry *
cell_get_dark_matter_sorts(const struct cell *c, const int sid) {
    
    /* We need to find at what position in the meta-array of
     sorts where the corresponding sid has been allocated since
     there might be gaps as we only allocated the directions that
     are in use.
     We create a mask with all the bits before the sid's one set to 1
     and apply it on the list of allocated directions. We then count
     the number of bits that are in the results to obtain the position
     of the correspondin sid in the meta-array */
    const int j = intrinsics_popcount(c->dark_matter.sort_allocated & ((1 << sid) - 1));
    
    /* Return the corresponding array */
    return &c->dark_matter.sort[j * (c->dark_matter.count + 1)];
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

  if (c->stars.sort != NULL) {
    swift_free("stars.sort", c->stars.sort);
    c->stars.sort = NULL;
    c->stars.sort_allocated = 0;
  }
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

/** Set the given flag for the given cell. */
__attribute__((always_inline)) INLINE static void cell_set_flag(struct cell *c,
                                                                uint32_t flag) {
  atomic_or(&c->flags, flag);
}

/** Clear the given flag for the given cell. */
__attribute__((always_inline)) INLINE static void cell_clear_flag(
    struct cell *c, uint32_t flag) {
  atomic_and(&c->flags, ~flag);
}

/** Get the given flag for the given cell. */
__attribute__((always_inline)) INLINE static int cell_get_flag(
    const struct cell *c, uint32_t flag) {
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

#endif /* SWIFT_CELL_H */
