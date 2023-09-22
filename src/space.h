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
#ifndef SWIFT_SPACE_H
#define SWIFT_SPACE_H

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <stddef.h>

/* Includes. */
#include "hydro_space.h"
#include "lock.h"
#include "parser.h"
#include "part.h"
#include "space_unique_id.h"
#include "velociraptor_struct.h"

/* Avoid cyclic inclusions */
struct cell;
struct cosmology;
struct gravity_props;
struct star_formation;
struct hydro_props;

/* Some constants. */
#define space_cellallocchunk 1000
#define space_splitsize_default 400
#define space_maxsize_default 8000000
#define space_extra_parts_default 0
#define space_extra_gparts_default 0
#define space_extra_sparts_default 100
#define space_extra_bparts_default 0
#define space_extra_sinks_default 0
#define space_expected_max_nr_strays_default 100
#define space_subsize_pair_hydro_default 256000000
#define space_subsize_self_hydro_default 32000
#define space_subsize_pair_stars_default 256000000
#define space_subsize_self_stars_default 32000
#define space_subsize_pair_grav_default 256000000
#define space_subsize_self_grav_default 32000
#define space_subdepth_diff_grav_default 4
#define space_max_top_level_cells_default 12
#define space_stretch 1.10f
#define space_maxreldx 0.1f

/* Maximum allowed depth of cell splits. */
#define space_cell_maxdepth 52

/* Globals needed in contexts without a space struct. Remember to dump and
 * restore these. */
extern int space_splitsize;
extern int space_maxsize;
extern int space_subsize_pair_hydro;
extern int space_subsize_self_hydro;
extern int space_subsize_pair_stars;
extern int space_subsize_self_stars;
extern int space_subsize_pair_grav;
extern int space_subsize_self_grav;
extern int space_subdepth_diff_grav;
extern int space_extra_parts;
extern int space_extra_gparts;
extern int space_extra_sparts;
extern int space_extra_bparts;
extern int space_extra_sinks;
extern double engine_redistribute_alloc_margin;
extern double engine_foreign_alloc_margin;

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

  /*! Are we doing star formation through sink particles? */
  int with_sink;

  /*! Are we running with some regular DM particles? */
  int with_DM;

  /*! Are we running with some DM background particles? */
  int with_DM_background;

  /*! Are we running with neutrino particles? */
  int with_neutrinos;

  /*! Width of the top-level cells. */
  double width[3];

  /*! Inverse of the top-level cell width */
  double iwidth[3];

  /*! The minimum top-level cell width allowed. */
  double cell_min;

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

  /*! Buffer of unused cells for the sub-cells. One chunk per thread. */
  struct cell **cells_sub;

  /*! The multipoles associated with the top-level (level 0) cells */
  struct gravity_tensors *multipoles_top;

  /*! Buffer of unused multipoles for the sub-cells. One chunk per thread. */
  struct gravity_tensors **multipoles_sub;

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

  /*! The total number of #bpart in the space. */
  size_t nr_bparts;

  /*! The total number of #sink in the space. */
  size_t nr_sinks;

  /*! The total number of neutrino #gpart in the space. */
  size_t nr_nuparts;

  /*! The total number of #part we allocated memory for */
  size_t size_parts;

  /*! The total number of #gpart we allocated memory for */
  size_t size_gparts;

  /*! The total number of #spart we allocated memory for */
  size_t size_sparts;

  /*! The total number of #bpart we allocated memory for */
  size_t size_bparts;

  /*! The total number of #sink we allocated memory for. */
  size_t size_sinks;

  /*! Number of inhibted gas particles in the space */
  size_t nr_inhibited_parts;

  /*! Number of inhibted gravity particles in the space */
  size_t nr_inhibited_gparts;

  /*! Number of inhibted star particles in the space */
  size_t nr_inhibited_sparts;

  /*! Number of inhibted black hole particles in the space */
  size_t nr_inhibited_bparts;

  /*! Number of inhibted sinks in the space */
  size_t nr_inhibited_sinks;

  /*! Number of extra #part we allocated (for on-the-fly creation) */
  size_t nr_extra_parts;

  /*! Number of extra #gpart we allocated (for on-the-fly creation) */
  size_t nr_extra_gparts;

  /*! Number of extra #spart we allocated (for on-the-fly creation) */
  size_t nr_extra_sparts;

  /*! Number of extra #bpart we allocated (for on-the-fly creation) */
  size_t nr_extra_bparts;

  /*! Number of extra #sink we allocated (for on-the-fly creation) */
  size_t nr_extra_sinks;

  /*! The particle data (cells have pointers to this). */
  struct part *parts;

  /*! The extended particle data (cells have pointers to this). */
  struct xpart *xparts;

  /*! The g-particle data (cells have pointers to this). */
  struct gpart *gparts;

  /*! The s-particle data (cells have pointers to this). */
  struct spart *sparts;

  /*! The b-particle data (cells have pointers to this). */
  struct bpart *bparts;

  /*! The sink particle data (cells have pointers to this). */
  struct sink *sinks;

  /*! Minimal mass of all the #part */
  float min_part_mass;

  /*! Minimal mass of all the dark-matter #gpart */
  float min_gpart_mass;

  /*! Minimal mass of all the #spart */
  float min_spart_mass;

  /*! Minimal mass of all the #sink */
  float min_sink_mass;

  /*! Minimal mass of all the #bpart */
  float min_bpart_mass;

  /*! Sum of the norm of the velocity of all the #part */
  float sum_part_vel_norm;

  /*! Sum of the norm of the velocity of all the dark-matter #gpart */
  float sum_gpart_vel_norm;

  /*! Sum of the norm of the velocity of all the #spart */
  float sum_spart_vel_norm;

  /*! Sum of the norm of the velocity of all the #sink */
  float sum_sink_vel_norm;

  /*! Sum of the norm of the velocity of all the #bpart */
  float sum_bpart_vel_norm;

  /*! Minimal gravity acceleration accross all particles */
  float min_a_grav;

  /*! Max gravity softening accross all particles */
  float max_softening;

  /*! Max multipole power accross all top-level cells */
  float max_mpole_power[SELF_GRAVITY_MULTIPOLE_ORDER + 1];

  /* Initial mean mass of each particle type in the system. */
  double initial_mean_mass_particles[swift_type_count];

  /* Initial count of each particle type in the system. */
  long long initial_count_particles[swift_type_count];

  /*! Initial shift that was applied to all particles upon start-up. */
  double initial_shift[3];

  /*! Initial value of the smoothing length read from the parameter file */
  float initial_spart_h;

  /*! Initial value of the smoothing length read from the parameter file */
  float initial_bpart_h;

  /*! General-purpose lock for this space. */
  swift_lock_type lock;

  /*! Number of queues in the system. */
  int nr_queues;

  /*! The associated engine. */
  struct engine *e;

  /*! The group information returned by VELOCIraptor for each #gpart. */
  struct velociraptor_gpart_data *gpart_group_data;

  /*! Structure dealing with the computation of a unique ID */
  struct unique_id unique_id;

#ifdef WITH_MPI

  /*! Buffers for parts that we will receive from foreign cells. */
  struct part *parts_foreign;
  size_t nr_parts_foreign, size_parts_foreign;

  /*! Buffers for g-parts that we will receive from foreign cells. */
  struct gpart *gparts_foreign;
  size_t nr_gparts_foreign, size_gparts_foreign;

  /*! Buffers for s-parts that we will receive from foreign cells. */
  struct spart *sparts_foreign;
  size_t nr_sparts_foreign, size_sparts_foreign;

  /*! Buffers for b-parts that we will receive from foreign cells. */
  struct bpart *bparts_foreign;
  size_t nr_bparts_foreign, size_bparts_foreign;

#endif
};

/* Function prototypes. */
void space_free_buff_sort_indices(struct space *s);
void space_parts_sort(struct part *parts, struct xpart *xparts, int *ind,
                      int *counts, int num_bins, ptrdiff_t parts_offset);
void space_gparts_sort(struct gpart *gparts, struct part *parts,
                       struct sink *sinks, struct spart *sparts,
                       struct bpart *bparts, int *ind, int *counts,
                       int num_bins);
void space_sparts_sort(struct spart *sparts, int *ind, int *counts,
                       int num_bins, ptrdiff_t sparts_offset);
void space_bparts_sort(struct bpart *bparts, int *ind, int *counts,
                       int num_bins, ptrdiff_t bparts_offset);
void space_sinks_sort(struct sink *sinks, int *ind, int *counts, int num_bins,
                      ptrdiff_t sinks_offset);
void space_getcells(struct space *s, int nr_cells, struct cell **cells,
                    const short int tid);
void space_init(struct space *s, struct swift_params *params,
                const struct cosmology *cosmo, double dim[3],
                const struct hydro_props *hydro_properties, struct part *parts,
                struct gpart *gparts, struct sink *sinks, struct spart *sparts,
                struct bpart *bparts, size_t Npart, size_t Ngpart, size_t Nsink,
                size_t Nspart, size_t Nbpart, size_t Nnupart, int periodic,
                int replicate, int remap_ids, int generate_gas_in_ics,
                int hydro, int gravity, int star_formation, int with_sink,
                int with_DM, int with_DM_background, int neutrinos, int verbose,
                int dry_run, int nr_nodes);
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
void space_regrid(struct space *s, int verbose);
void space_allocate_extras(struct space *s, int verbose);
void space_split(struct space *s, int verbose);
void space_reorder_extras(struct space *s, int verbose);
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
void space_bparts_get_cell_index(struct space *s, int *sind, int *cell_counts,
                                 size_t *count_inhibited_bparts,
                                 size_t *count_extra_bparts, int verbose);
void space_sinks_get_cell_index(struct space *s, int *sind, int *cell_counts,
                                size_t *count_inhibited_sinks,
                                size_t *count_extra_sinks, int verbose);
void space_synchronize_particle_positions(struct space *s);
void space_first_init_parts(struct space *s, int verbose);
void space_first_init_gparts(struct space *s, int verbose);
void space_first_init_sparts(struct space *s, int verbose);
void space_first_init_bparts(struct space *s, int verbose);
void space_first_init_sinks(struct space *s, int verbose);
void space_collect_mean_masses(struct space *s, int verbose);
void space_init_parts(struct space *s, int verbose);
void space_init_gparts(struct space *s, int verbose);
void space_init_sparts(struct space *s, int verbose);
void space_init_bparts(struct space *s, int verbose);
void space_init_sinks(struct space *s, int verbose);
void space_after_snap_tracer(struct space *s, int verbose);
void space_convert_quantities(struct space *s, int verbose);
void space_convert_rt_quantities(struct space *s, int verbose);
void space_post_init_parts(struct space *s, int verbose);
void space_link_cleanup(struct space *s);
void space_check_drift_point(struct space *s, integertime_t ti_drift,
                             int multipole);
void space_check_top_multipoles_drift_point(struct space *s,
                                            integertime_t ti_drift);
void space_check_timesteps(const struct space *s);
void space_check_limiter(struct space *s);
void space_check_swallow(struct space *s);
void space_check_sort_flags(struct space *s);
void space_remap_ids(struct space *s, int nr_nodes, int verbose);
long long space_get_max_parts_id(struct space *s);
void space_replicate(struct space *s, int replicate, int verbose);
void space_generate_gas(struct space *s, const struct cosmology *cosmo,
                        const struct hydro_props *hydro_properties,
                        const int periodic, const int with_DM_background,
                        const int with_neutrinos, const double dim[3],
                        const int verbose);
void space_check_cosmology(struct space *s, const struct cosmology *cosmo,
                           const int with_hydro, const int rank,
                           const int check_neutrinos);
void space_reset_task_counters(struct space *s);
void space_clean(struct space *s);
void space_free_cells(struct space *s);

void space_free_foreign_parts(struct space *s, const int clear_cell_pointers);

void space_struct_dump(struct space *s, FILE *stream);
void space_struct_restore(struct space *s, FILE *stream);
void space_write_cell_hierarchy(const struct space *s, int j);
void space_compute_star_formation_stats(const struct space *s,
                                        struct star_formation *star_form);
void space_check_unskip_flags(const struct space *s);
#endif /* SWIFT_SPACE_H */
