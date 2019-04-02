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
#ifndef SWIFT_SPACE_H
#define SWIFT_SPACE_H

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stddef.h>

/* Includes. */
#include "gravity_properties.h"
#include "hydro_space.h"
#include "lock.h"
#include "parser.h"
#include "part.h"
#include "velociraptor_struct.h"

/* Avoid cyclic inclusions */
struct cell;
struct cosmology;

/* Some constants. */
#define space_cellallocchunk 1000
#define space_splitsize_default 400
#define space_maxsize_default 8000000
#define space_extra_parts_default 0
#define space_extra_gparts_default 0
#define space_extra_sparts_default 100
#define space_expected_max_nr_strays_default 100
#define space_subsize_pair_hydro_default 256000000
#define space_subsize_self_hydro_default 32000
#define space_subsize_pair_grav_default 256000000
#define space_subsize_self_grav_default 32000
#define space_subdepth_diff_grav_default 4
#define space_max_top_level_cells_default 12
#define space_stretch 1.10f
#define space_maxreldx 0.1f

/* Maximum allowed depth of cell splits. */
#define space_cell_maxdepth 52

/* Split size. */
extern int space_splitsize;
extern int space_maxsize;
extern int space_subsize_pair_hydro;
extern int space_subsize_self_hydro;
extern int space_subsize_pair_grav;
extern int space_subsize_self_grav;
extern int space_subsize_pair_stars;
extern int space_subsize_self_stars;
extern int space_subdepth_diff_grav;
extern int space_extra_parts;
extern int space_extra_gparts;
extern int space_extra_sparts;

/**
 * @brief The space in which the cells and particles reside.
 */
struct space {

  /*! Spatial extent. */
  double dim[3];

  /*! Is the space periodic? */
  int periodic;

  /*! Extra space information needed for some hydro schemes. */
  struct hydro_space hs;

  /*! Are we doing hydrodynamics? */
  int with_hydro;

  /*! Are we doing gravity? */
  int with_self_gravity;

  /*! Are we doing star formation? */
  int with_star_formation;

  /*! Width of the top-level cells. */
  double width[3];

  /*! Inverse of the top-level cell width */
  double iwidth[3];

  /*! The minimum top-level cell width allowed. */
  double cell_min;

  /*! Current maximum displacement for particles. */
  float dx_max;

  /*! Space dimensions in number of top-cells. */
  int cdim[3];

  /*! Maximal depth reached by the tree */
  int maxdepth;

  /*! Number of top-level cells. */
  int nr_cells;

  /*! Total number of cells (top- and sub-) */
  int tot_cells;

  /*! Number of *local* top-level cells */
  int nr_local_cells;

  /*! Number of *local* top-level cells with tasks */
  int nr_local_cells_with_tasks;

  /*! Number of top-level cells that have >0 particle (of any kind) */
  int nr_cells_with_particles;

  /*! Number of top-level cells that have >0 particle (of any kind) */
  int nr_local_cells_with_particles;

  /*! The (level 0) cells themselves. */
  struct cell *cells_top;

  /*! Buffer of unused cells for the sub-cells. */
  struct cell *cells_sub;

  /*! The multipoles associated with the top-level (level 0) cells */
  struct gravity_tensors *multipoles_top;

  /*! Buffer of unused multipoles for the sub-cells. */
  struct gravity_tensors *multipoles_sub;

  /*! The indices of the *local* top-level cells */
  int *local_cells_top;

  /*! The indices of the *local* top-level cells with tasks */
  int *local_cells_with_tasks_top;

  /*! The indices of the top-level cells that have >0 particles (of any kind) */
  int *cells_with_particles_top;

  /*! The indices of the top-level cells that have >0 particles (of any kind) */
  int *local_cells_with_particles_top;

  /*! The total number of #part in the space. */
  size_t nr_parts;

  /*! The total number of #gpart in the space. */
  size_t nr_gparts;

  /*! The total number of #spart in the space. */
  size_t nr_sparts;

  /*! The total number of #part we allocated memory for */
  size_t size_parts;

  /*! The total number of #gpart we allocated memory for */
  size_t size_gparts;

  /*! The total number of #spart we allocated memory for */
  size_t size_sparts;

  /*! Number of inhibted gas particles in the space */
  size_t nr_inhibited_parts;

  /*! Number of inhibted gravity particles in the space */
  size_t nr_inhibited_gparts;

  /*! Number of inhibted star particles in the space */
  size_t nr_inhibited_sparts;

  /*! Number of extra #part we allocated (for on-the-fly creation) */
  size_t nr_extra_parts;

  /*! Number of extra #gpart we allocated (for on-the-fly creation) */
  size_t nr_extra_gparts;

  /*! Number of extra #spart we allocated (for on-the-fly creation) */
  size_t nr_extra_sparts;

  /*! The particle data (cells have pointers to this). */
  struct part *parts;

  /*! The extended particle data (cells have pointers to this). */
  struct xpart *xparts;

  /*! The g-particle data (cells have pointers to this). */
  struct gpart *gparts;

  /*! The s-particle data (cells have pointers to this). */
  struct spart *sparts;

  /*! Minimal mass of all the #part */
  float min_part_mass;

  /*! Minimal mass of all the dark-matter #gpart */
  float min_gpart_mass;

  /*! Minimal mass of all the #spart */
  float min_spart_mass;

  /*! Sum of the norm of the velocity of all the #part */
  float sum_part_vel_norm;

  /*! Sum of the norm of the velocity of all the dark-matter #gpart */
  float sum_gpart_vel_norm;

  /*! Sum of the norm of the velocity of all the #spart */
  float sum_spart_vel_norm;

  /*! Initial value of the smoothing length read from the parameter file */
  float initial_spart_h;

  /*! General-purpose lock for this space. */
  swift_lock_type lock;

  /*! Number of queues in the system. */
  int nr_queues;

  /*! The associated engine. */
  struct engine *e;

  /*! The FOF group data. */
  struct fof fof_data;

  /*! The group information returned by VELOCIraptor for each #gpart. */
  struct velociraptor_gpart_data *gpart_group_data;

#ifdef WITH_MPI

  /*! Buffers for parts that we will receive from foreign cells. */
  struct part *parts_foreign;
  size_t nr_parts_foreign, size_parts_foreign;

  /*! Buffers for g-parts that we will receive from foreign cells. */
  struct gpart *gparts_foreign;
  size_t nr_gparts_foreign, size_gparts_foreign;

  /*! Buffers for g-parts that we will receive from foreign cells. */
  struct spart *sparts_foreign;
  size_t nr_sparts_foreign, size_sparts_foreign;

#endif
};

/* Function prototypes. */
void space_free_buff_sort_indices(struct space *s);
void space_parts_sort(struct part *parts, struct xpart *xparts, int *ind,
                      int *counts, int num_bins, ptrdiff_t parts_offset);
void space_gparts_sort(struct gpart *gparts, struct part *parts,
                       struct spart *sparts, int *ind, int *counts,
                       int num_bins);
void space_sparts_sort(struct spart *sparts, int *ind, int *counts,
                       int num_bins, ptrdiff_t sparts_offset);
void space_getcells(struct space *s, int nr_cells, struct cell **cells);
void space_init(struct space *s, struct swift_params *params,
                const struct cosmology *cosmo, double dim[3],
                struct part *parts, struct gpart *gparts, struct spart *sparts,
                size_t Npart, size_t Ngpart, size_t Nspart, int periodic,
                int replicate, int generate_gas_in_ics, int hydro, int gravity,
                int star_formation, int verbose, int dry_run);
void space_sanitize(struct space *s);
void space_map_cells_pre(struct space *s, int full,
                         void (*fun)(struct cell *c, void *data), void *data);
void space_map_parts(struct space *s,
                     void (*fun)(struct part *p, struct cell *c, void *data),
                     void *data);
void space_map_parts_xparts(struct space *s,
                            void (*fun)(struct part *p, struct xpart *xp,
                                        struct cell *c));
void space_map_cells_post(struct space *s, int full,
                          void (*fun)(struct cell *c, void *data), void *data);
void space_rebuild(struct space *s, int repartitioned, int verbose);
void space_recycle(struct space *s, struct cell *c);
void space_recycle_list(struct space *s, struct cell *cell_list_begin,
                        struct cell *cell_list_end,
                        struct gravity_tensors *multipole_list_begin,
                        struct gravity_tensors *multipole_list_end);
void space_split(struct space *s, int verbose);
void space_reorder_extras(struct space *s, int verbose);
void space_split_mapper(void *map_data, int num_elements, void *extra_data);
void space_list_useful_top_level_cells(struct space *s);
void space_parts_get_cell_index(struct space *s, int *ind, int *cell_counts,
                                size_t *count_inhibited_parts,
                                size_t *count_extra_parts, int verbose);
void space_gparts_get_cell_index(struct space *s, int *gind, int *cell_counts,
                                 size_t *count_inhibited_gparts,
                                 size_t *count_extra_gparts, int verbose);
void space_sparts_get_cell_index(struct space *s, int *sind, int *cell_counts,
                                 size_t *count_inhibited_sparts,
                                 size_t *count_extra_sparts, int verbose);
void space_synchronize_particle_positions(struct space *s);
void space_do_parts_sort(void);
void space_do_gparts_sort(void);
void space_do_sparts_sort(void);
void space_first_init_parts(struct space *s, int verbose);
void space_first_init_gparts(struct space *s, int verbose);
void space_first_init_sparts(struct space *s, int verbose);
void space_init_parts(struct space *s, int verbose);
void space_init_gparts(struct space *s, int verbose);
void space_init_sparts(struct space *s, int verbose);
void space_convert_quantities(struct space *s, int verbose);
void space_link_cleanup(struct space *s);
void space_check_drift_point(struct space *s, integertime_t ti_drift,
                             int multipole);
void space_check_top_multipoles_drift_point(struct space *s,
                                            integertime_t ti_drift);
void space_check_timesteps(struct space *s);
void space_check_limiter(struct space *s);
void space_replicate(struct space *s, int replicate, int verbose);
void space_generate_gas(struct space *s, const struct cosmology *cosmo,
                        int periodic, const double dim[3], int verbose);
void space_check_cosmology(struct space *s, const struct cosmology *cosmo,
                           int rank);
void space_reset_task_counters(struct space *s);
void space_clean(struct space *s);
void space_free_cells(struct space *s);

void space_free_foreign_parts(struct space *s);

void space_struct_dump(struct space *s, FILE *stream);
void space_struct_restore(struct space *s, FILE *stream);
void space_write_cell_hierarchy(const struct space *s);

#endif /* SWIFT_SPACE_H */
