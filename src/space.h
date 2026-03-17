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
#include <fftw3.h>

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
#define space_grid_split_threshold_default 400
#define space_extra_parts_default 0
#define space_extra_gparts_default 200
#define space_extra_sparts_default 200
#define space_extra_bparts_default 0
#define space_extra_sinks_default 0
#define space_expected_max_nr_strays_default 100
#define space_subsize_pair_hydro_default 256000000
#define space_subsize_self_hydro_default 32000
#define space_subsize_pair_stars_default 256000000
#define space_subsize_self_stars_default 32000
#define space_subsize_pair_grav_default 256000000
#define space_subsize_self_grav_default 32000
#define space_recurse_size_self_hydro_default 100
#define space_recurse_size_pair_hydro_default 100
#define space_recurse_size_self_stars_default 100
#define space_recurse_size_pair_stars_default 100
#define space_recurse_size_self_black_holes_default 100
#define space_recurse_size_pair_black_holes_default 100
#define space_recurse_size_self_sinks_default 100
#define space_recurse_size_pair_sinks_default 100
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
extern int space_grid_split_threshold;
extern int space_subsize_pair_hydro;
extern int space_subsize_self_hydro;
extern int space_subsize_pair_stars;
extern int space_subsize_self_stars;
extern int space_subsize_pair_grav;
extern int space_subsize_self_grav;
extern int space_subdepth_diff_grav;
extern int space_recurse_size_self_hydro;
extern int space_recurse_size_pair_hydro;
extern int space_recurse_size_self_stars;
extern int space_recurse_size_pair_stars;
extern int space_recurse_size_self_black_holes;
extern int space_recurse_size_pair_black_holes;
extern int space_recurse_size_self_sinks;
extern int space_recurse_size_pair_sinks;
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

  /*! Have the top-level cells' time-steps been updated? */
  char *cells_top_updated;

#ifdef WITH_MPI

  /*! Buffers for parts that we will receive from foreign cells. */
  struct part *parts_foreign;
  size_t nr_parts_foreign, size_parts_foreign;

  /*! Buffers for g-parts that we will receive from foreign cells. */
  struct gpart_foreign *gparts_foreign;
  struct gpart_fof_foreign *gparts_fof_foreign;
  size_t nr_gparts_foreign, size_gparts_foreign;

  /*! Buffers for s-parts that we will receive from foreign cells. */
  struct spart *sparts_foreign;
  size_t nr_sparts_foreign, size_sparts_foreign;

  /*! Buffers for b-parts that we will receive from foreign cells. */
  struct bpart *bparts_foreign;
  size_t nr_bparts_foreign, size_bparts_foreign;

  /*! Buffers for sink-parts that we will receive from foreign cells. */
  struct sink *sinks_foreign;
  size_t nr_sinks_foreign, size_sinks_foreign;

#endif
};

struct cic_mapper_data;
struct tree_data;
//struct cell_basic;
struct AMR_levels {
  struct cell **cells;
  double mean_density; 
  int cell_count; 
  int ghost_count;
  int depth; 
  int cdim;
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
                    const short int tid, int index, int ghost);
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
void space_mark_cell_as_updated(struct space *s, const struct cell *c);
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
void space_get_density(struct space *s, struct swift_params *params, struct engine *e, int use_multigrid);
void space_apply_FMG(struct space *s, struct engine *e);
void gpart_to_mesh_CIC_mapper(void* map_data, int num, void* extra);
void mesh_to_gpart_CIC_mapper(void* map_data, int num, void* extra);
void density_to_cells(struct cic_mapper_data* data, struct cell* top_cells, int nr_cells, int cdim[3]);
void get_residual(double *pot, double *dens, int cdim[3], double multiplier, double *residual, double delta);
void get_residual_array(double *pot, double *dens, int cdim[3], double multiplier, double *residual, double delta);
void perform_red_black_sweep(double *pot, double *dens, int cdim[3], double multiplier, double delta);
void set_initial_guess(double *potential_array, int cdim[3]);
double get_mean_density(double *dens, const int N);
void get_pm_potential(struct cic_mapper_data* data, const int N, const double box_size, struct threadpool* tp, int cdim[3]);
void mesh_apply_Green_function(struct threadpool* tp, fftw_complex* frho,
                               const int slice_offset, const int slice_width,
                               const int N, const double r_s,
                               const double box_size, const int deconvolve, const int discrete_symbol);
void restrict_problem(double *H_array, double *residual, int cdim[3], double delta);
void solve_coarser_problem(double *pot, double *residual, int cdim[3], double delta, const int N_stop, const int N_start);
void prolongate_problem(double *coarser_solution, double *pot, int cdim[3]);
void apply_GS(double *density, double *pot, int cdim[3], double mean_density, double box_size);
void apply_multigrid(double *density, double *pot, int cdim[3], double mean_density, double box_size, int N_stop, int N_start, int V_max);
void prolongate_solution(double *pot_coarse, double *pot_fine, const int N, const int N_double);
void get_accelerations(struct cic_mapper_data* data, double *pot, struct threadpool* tp, struct space *s,const int N, int step);
void space_get_AMR_density(struct space *s, struct engine *e, int level_check);
void check_progeny(struct space *s, struct cell *parent, int *length, int *level);
void construct_daughter(struct cell *parent, int i, int *length);
void destroy_daughters(struct space *s, struct cell *cells_top, int nr_cells);
void level_down(struct space *s, struct cell *parent, int *length);
void assign_densities(struct cell *cells_top, struct threadpool *tp, int nr_cells, const double boxsize, int level_check);
void initialise_tree_search(struct gpart *part, struct cell *home_cell, struct cell *cells_top, const double boxsize, int nr_topcells, int level_check, int verbose);
void check_lower_level(struct cell *parent, struct cell *cells_top, const double boxsize, int nr_topcells, int *level, int level_check);
double get_overlap(double pbox[6], double cloc[3], double cwidth[3]);
void search_tree(struct gpart *part, double pbox[6], struct cell *cells_top, double width, int nr_topcells, double *acc, int level_check, int reverse, int verbose);
void assign_lower_level(struct gpart *part, double pbox[6], struct cell *parent, double width, int *level, double *acc, int level_check, int reverse, int verbose);
void assign_parent_densities(struct cell *parent, int *level);
void export_grid_data(struct cell *cells_top, int nr_cells);
void test_density_assignment(struct cell *cells_top, int nr_cells, const double boxsize, int nr_gparts, struct space *s, int level_check);
//void assign_CIC_densities(struct space *s);
//void assign_densities_test(struct space *s);
//void search_tree_test(struct gpart *part, double pbox[6], struct cell *cells_top, double width, int nr_topcells, double *CIC_overlap, int *CIC_cell, int *counter);
void get_CIC_results(struct gpart *gp, struct cell *cells_top, struct space *s);
//void get_AMR_results(struct gpart *gp, struct cell *cells_top, struct space *s, double *CIC_overlap, int *CIC_cell);
int perform_uniform_calculation(struct space *s, int min_depth, int max_depth, struct AMR_levels levels[max_depth+1]);
//int perform_uniform_calculation(struct space *s, struct cell *cells_top, int N_levels);
void sort_lower_level(struct cell *parent, double *rho, double fac, int max_level, int *current_level, int cdim[3]);
void potential_to_cells(struct cell *cells_top, double *pot, int grid_size, double grid_top, double fac, int level);
void to_lower_level(struct cell *cell, double pot_cell[6], double pot, int level, int *current_level);
void cells_top_mapper(void* map_data, const int num, void* extra);
//void perform_nonuniform_calculation(struct space *s, int min_depth, int max_depth, struct AMR_levels levels[max_depth+1], int max_gridsize);
void perform_nonuniform_calculation(struct space *s, int min_depth, int max_depth, struct AMR_levels *levels, int max_gridsize, double box_size);
void extract_AMR_patches(struct space *s, struct cell *cells_top, int extract_depth, struct cell ***extracted_cells, int *count);
void extract_lower(struct cell *parent, int *current_depth, int extract_depth, struct cell ***extracted_cells, int *count);
int sort_cell(const void *a, const void *b);
void set_patch_guess(struct space *s, struct AMR_levels *coarse, struct AMR_levels *fine, int nr_cells, int gridsize, double box_size, int min_depth); 
double find_coarse_cell(struct AMR_levels *coarse, int cid, int missing, double fac_coarse, int cdimH[3]);
int find_neighbour(struct AMR_levels *fine, struct cell *parent, struct cell *cell_link, int search_id, double fac, int which_neighbour, int pre_smoothing);
//void create_ghost(struct space *s, struct AMR_levels *coarse, struct AMR_levels *fine, double search_loc[3], int which_neighbour, double fac_coarse);
void create_ghost(struct space *s, struct AMR_levels *coarse, struct AMR_levels *fine, struct cell *parent, size_t search_id, double fac_coarse, int min_depth, int which_neighbour);
void initialise_link_neighbours(struct AMR_levels *fine, double box_size, int gridsize);
void link_neighbour(struct cell *link_cell, struct AMR_levels *fine, double search_loc[3], double fac, int which_neighbour);
void extrapolate_mask_values(struct AMR_levels *coarse);
void get_patch_density(struct space *s, struct AMR_levels *level);
void perform_patch_sweep(struct AMR_levels *level, double delta);
void transfer_residual_array(struct AMR_levels fine, struct AMR_levels *coarse, double delta);
void to_coarser_patch(struct AMR_levels *levels, double delta, int current_depth, int min_depth);
void perform_multigrid_acceleration(struct space *s, int min_depth, int max_depth, struct AMR_levels levels[max_depth+1], int current_depth);
double get_patch_residual(struct AMR_levels level, double delta);
//void replace_masked_cells(struct AMR_levels *level, double delta);
void perform_masked_patch_sweep(struct AMR_levels *level, double delta);
double get_masked_residual(struct AMR_levels level, double delta, int nr_active, int step_nr);
void link_uniform_level(struct AMR_levels *level);
double get_solution_magnitude(struct AMR_levels level);
void converge_first_order(struct AMR_levels *level, double delta, int nr_steps);
void first_order_sweep(struct AMR_levels *level, double delta);
double get_first_order_residual(struct AMR_levels *level, double delta);
void to_finer_patch(struct AMR_levels *levels, int target_depth);
void transfer_coarser_residual(struct AMR_levels fine, struct AMR_levels *coarse, double delta);
void build_refinement_map(struct space *s, int min_depth, int max_depth, struct AMR_levels levels[max_depth+1]);
void modify_tree(struct space *s, int min_depth, int max_depth, struct AMR_levels levels[max_depth+1], int new_cell_count[max_depth +1]);
void divide_particles(struct cell *c, struct space *s);
//void create_concept_link(struct space *s, struct AMR_levels *level, struct cell **cells, int nr_cells);
int search_neighbours(struct AMR_levels *level, double search_loc[3], double fac);
void mark_neighbours(struct space *s, int min_depth, struct AMR_levels *level, struct cell *curr_cell);
void create_link(struct space *s, struct AMR_levels *level, struct cell **cells, int nr_cells, int which_neighbour);
void potential_to_gparts(struct space *s, int min_depth, int max_depth, struct AMR_levels *levels);
void get_AMR_potential(struct space *s, int max_depth, int current_depth, struct AMR_levels *levels, struct gpart *gp, int cell_nr);
double CIC_get_AMR(struct space *s, struct gpart *gp, double x[3], double width[3], double boxsize);
void link_nonuniform_level(struct space *s, struct AMR_levels *level, int start_index, int link_nr);
void interpolate_trilinear(struct AMR_levels *coarse, struct AMR_levels *fine);
void free_gparts_in_cells(struct cell *c, int *level);
void init_test_single_particle(struct engine *e);
void get_progeny(struct space *s, struct cell *c, int deisred_depth, int *curr_depth);
void particle_to_cells_recursive(double part_loc[3], struct cell **cells, int nr_cells);
void generate_particles(struct space *s, int N_parts_new, double positions[N_parts_new][3], double r_parts);
void potential_to_fake_gparts(struct space *s, int min_depth, int max_depth, struct AMR_levels levels[max_depth+1]);
#endif /* SWIFT_SPACE_H */
