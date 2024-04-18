/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *                    Angus Lepper (angus.lepper@ed.ac.uk)
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
#ifndef SWIFT_ENGINE_H
#define SWIFT_ENGINE_H

/* Config parameters. */
#include <config.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Includes. */
#include "barrier.h"
#include "clocks.h"
#include "collectgroup.h"
#include "ic_info.h"
#include "lightcone/lightcone.h"
#include "lightcone/lightcone_array.h"
#include "mesh_gravity.h"
#include "output_options.h"
#include "parser.h"
#include "partition.h"
#include "runner.h"
#include "scheduler.h"
#include "space.h"
#include "task.h"
#include "tracers_triggers.h"
#include "units.h"
#include "velociraptor_interface.h"

struct black_holes_properties;
struct extra_io_properties;
struct external_potential;
struct forcing_terms;

/**
 * @brief The different policies the #engine can follow.
 */
enum engine_policy {
  engine_policy_none = 0,
  engine_policy_rand = (1 << 0),
  engine_policy_steal = (1 << 1),
  engine_policy_keep = (1 << 2),
  engine_policy_block = (1 << 3),
  engine_policy_cputight = (1 << 4),
  engine_policy_mpi = (1 << 5),
  engine_policy_setaffinity = (1 << 6),
  engine_policy_hydro = (1 << 7),
  engine_policy_self_gravity = (1 << 8),
  engine_policy_external_gravity = (1 << 9),
  engine_policy_cosmology = (1 << 10),
  engine_policy_drift_all = (1 << 11),
  engine_policy_reconstruct_mpoles = (1 << 12),
  engine_policy_temperature = (1 << 13),
  engine_policy_cooling = (1 << 14),
  engine_policy_stars = (1 << 15),
  engine_policy_structure_finding = (1 << 16),
  engine_policy_star_formation = (1 << 17),
  engine_policy_feedback = (1 << 18),
  engine_policy_black_holes = (1 << 19),
  engine_policy_fof = (1 << 20),
  engine_policy_timestep_limiter = (1 << 21),
  engine_policy_timestep_sync = (1 << 22),
  engine_policy_csds = (1 << 23),
  engine_policy_line_of_sight = (1 << 24),
  engine_policy_sinks = (1 << 25),
  engine_policy_rt = (1 << 26),
  engine_policy_power_spectra = (1 << 27),
};
#define engine_maxpolicy 28
extern const char *engine_policy_names[engine_maxpolicy + 1];

/**
 * @brief The different unusual events that can take place in a time-step.
 */
enum engine_step_properties {
  engine_step_prop_none = 0,
  engine_step_prop_rebuild = (1 << 0),
  engine_step_prop_redistribute = (1 << 1),
  engine_step_prop_repartition = (1 << 2),
  engine_step_prop_statistics = (1 << 3),
  engine_step_prop_snapshot = (1 << 4),
  engine_step_prop_restarts = (1 << 5),
  engine_step_prop_stf = (1 << 6),
  engine_step_prop_fof = (1 << 7),
  engine_step_prop_mesh = (1 << 8),
  engine_step_prop_power_spectra = (1 << 9),
  engine_step_prop_done = (1 << 10),
};

/* Some constants */
#define engine_maxproxies 64
#define engine_tasksreweight 1
#define engine_parts_size_grow 1.05
#define engine_redistribute_alloc_margin_default 1.2
#define engine_rebuild_link_alloc_margin 1.2
#define engine_foreign_alloc_margin_default 1.05
#define engine_default_energy_file_name "statistics"
#define engine_default_timesteps_file_name "timesteps"
#define engine_default_rt_subcycles_file_name "rtsubcycles"
#define engine_max_parts_per_ghost_default 1000
#define engine_max_sparts_per_ghost_default 1000
#define engine_max_parts_per_cooling_default 10000
#define engine_star_resort_task_depth_default 2
#define engine_tasks_per_cell_margin 1.2
#define engine_default_stf_subdir_per_output "."
#define engine_default_snapshot_subdir "."

/**
 * @brief The rank of the engine as a global variable (for messages).
 */
extern int engine_rank;

/**
 * @brief The current step as a global variable (for messages).
 */
extern int engine_current_step;

/* Data structure for the engine. */
struct engine {

  /* Number of task threads on which to run. */
  int nr_threads;

  /* Number of threadpool threads on which to run. */
  int nr_pool_threads;

  /* The space with which the runner is associated. */
  struct space *s;

  /* The runner's threads. */
  struct runner *runners;

  /* The running policy. */
  int policy;

  /* The task scheduler. */
  struct scheduler sched;

  /* Common threadpool for all the engine's tasks. */
  struct threadpool threadpool;

  /* The minimum and maximum allowed dt */
  double dt_min, dt_max;

  /* Maximum time-step allowed by the RMS condition in cosmology runs. */
  double dt_max_RMS_displacement;

  /* Dimensionless factor for the RMS time-step condition. */
  double max_RMS_displacement_factor;

  /* When computing the max RMS dt, should only the gas particles
   * be considered as the baryon component? */
  int max_RMS_dt_use_only_gas;

  /* Time of the simulation beginning */
  double time_begin;

  /* Time of the simulation end */
  double time_end;

  /* The previous system time. */
  double time_old;
  integertime_t ti_old;

  /* The current system time. */
  double time;
  integertime_t ti_current;

  /* The earliest time any particle may still need to be drifted from */
  integertime_t ti_earliest_undrifted;

  /* The highest active bin at this time */
  timebin_t max_active_bin;

  /* The lowest active bin at this time */
  timebin_t min_active_bin;

  /* RT sub-cycling counterparts for timestepping vars */
  integertime_t ti_current_subcycle;
  timebin_t max_active_bin_subcycle;
  timebin_t min_active_bin_subcycle;

  /* Maximal number of radiative transfer sub-cycles per hydro step */
  int max_nr_rt_subcycles;

  /* Time step */
  double time_step;

  /* Time base */
  double time_base;
  double time_base_inv;

  /* Minimal hydro ti_end for the next time-step */
  integertime_t ti_hydro_end_min;

  /* Maximal hydro ti_end for the next time-step */
  integertime_t ti_hydro_end_max;

  /* Maximal hydro ti_beg for the next time-step */
  integertime_t ti_hydro_beg_max;

  /* Minimal rt ti_end for the next time-step */
  integertime_t ti_rt_end_min;

  /* Maximal rt ti_beg for the next time-step */
  integertime_t ti_rt_beg_max;

  /* Minimal gravity ti_end for the next time-step */
  integertime_t ti_gravity_end_min;

  /* Maximal gravity ti_end for the next time-step */
  integertime_t ti_gravity_end_max;

  /* Maximal gravity ti_beg for the next time-step */
  integertime_t ti_gravity_beg_max;

  /* Minimal stars ti_end for the next time-step */
  integertime_t ti_stars_end_min;

  /* Maximal stars ti_end for the next time-step */
  integertime_t ti_stars_end_max;

  /* Maximal stars ti_beg for the next time-step */
  integertime_t ti_stars_beg_max;

  /* Minimal black holes ti_end for the next time-step */
  integertime_t ti_black_holes_end_min;

  /* Maximal black holes ti_end for the next time-step */
  integertime_t ti_black_holes_end_max;

  /* Maximal black holes ti_beg for the next time-step */
  integertime_t ti_black_holes_beg_max;

  /* Minimal sinks ti_end for the next time-step */
  integertime_t ti_sinks_end_min;

  /* Maximal sinks ti_end for the next time-step */
  integertime_t ti_sinks_end_max;

  /* Maximal sinks ti_beg for the next time-step */
  integertime_t ti_sinks_beg_max;

  /* Minimal overall ti_end for the next time-step */
  integertime_t ti_end_min;

  /* Maximal overall ti_beg for the next time-step */
  integertime_t ti_beg_max;

  /* Number of particles updated in the previous step */
  long long updates, g_updates, s_updates, b_updates, sink_updates, rt_updates;

  /* Number of updates since the last rebuild */
  long long updates_since_rebuild;
  long long g_updates_since_rebuild;
  long long s_updates_since_rebuild;
  long long sink_updates_since_rebuild;
  long long b_updates_since_rebuild;

  /* Star formation logger information */
  struct star_formation_history_accumulator sfh;

  /* Properties of the previous step */
  int step_props;

  /* Total numbers of particles in the system. */
  long long total_nr_parts;
  long long total_nr_gparts;
  long long total_nr_sparts;
  long long total_nr_sinks;
  long long total_nr_bparts;
  long long total_nr_DM_background_gparts;
  long long total_nr_neutrino_gparts;

  /* Total numbers of cells (top-level and sub-cells) in the system. */
  long long total_nr_cells;

  /* Total numbers of tasks in the system. */
  long long total_nr_tasks;

  /* The total number of inhibited particles in the system. */
  long long nr_inhibited_parts;
  long long nr_inhibited_gparts;
  long long nr_inhibited_sparts;
  long long nr_inhibited_sinks;
  long long nr_inhibited_bparts;

#ifdef SWIFT_DEBUG_CHECKS
  /* Total number of particles removed from the system since the last rebuild */
  long long count_inhibited_parts;
  long long count_inhibited_gparts;
  long long count_inhibited_sparts;
  long long count_inhibited_bparts;
#endif

  /* Maximal ID of the parts, *excluding* background particles
   * (used for the generation of new IDs when splitting) */
  long long max_parts_id;

  /* Total mass in the simulation */
  double total_mass;

  /* Conversion factor between microscopic neutrino mass (eV) and gpart mass */
  double neutrino_mass_conversion_factor;

  /* The internal system of units */
  const struct unit_system *internal_units;

  /* Snapshot information */
  double a_first_snapshot;
  double time_first_snapshot;
  double delta_time_snapshot;

  /* Output_List for the snapshots */
  struct output_list *output_list_snapshots;

  /* Integer time of the next snapshot */
  integertime_t ti_next_snapshot;

  char snapshot_base_name[PARSER_MAX_LINE_SIZE];
  char snapshot_subdir[PARSER_MAX_LINE_SIZE];
  char snapshot_dump_command[PARSER_MAX_LINE_SIZE];
  int snapshot_subsample[swift_type_count];
  float snapshot_subsample_fraction[swift_type_count];
  int snapshot_run_on_dump;
  int snapshot_distributed;
  int snapshot_lustre_OST_count;
  int snapshot_compression;
  int snapshot_invoke_stf;
  int snapshot_invoke_fof;
  int snapshot_invoke_ps;
  struct unit_system *snapshot_units;
  int snapshot_use_delta_from_edge;
  double snapshot_delta_from_edge;
  int snapshot_output_count;

  /* Snapshot recording trigger mechanism counters */
  double snapshot_recording_triggers_part[num_snapshot_triggers_part];
  double snapshot_recording_triggers_desired_part[num_snapshot_triggers_part];
  int snapshot_recording_triggers_started_part[num_snapshot_triggers_part];

  double snapshot_recording_triggers_spart[num_snapshot_triggers_spart];
  double snapshot_recording_triggers_desired_spart[num_snapshot_triggers_spart];
  int snapshot_recording_triggers_started_spart[num_snapshot_triggers_spart];

  double snapshot_recording_triggers_bpart[num_snapshot_triggers_bpart];
  double snapshot_recording_triggers_desired_bpart[num_snapshot_triggers_bpart];
  int snapshot_recording_triggers_started_bpart[num_snapshot_triggers_bpart];

  /* Metadata from the ICs */
  struct ic_info *ics_metadata;

  /* Structure finding information */
  double a_first_stf_output;
  double time_first_stf_output;
  double delta_time_stf;

  /* Output_List for the structure finding */
  struct output_list *output_list_stf;

  /* Integer time of the next stf output */
  integertime_t ti_next_stf;

  char stf_config_file_name[PARSER_MAX_LINE_SIZE];
  char stf_base_name[PARSER_MAX_LINE_SIZE];
  char stf_subdir_per_output[PARSER_MAX_LINE_SIZE];
  int stf_output_count;

  /* FoF black holes seeding information */
  double a_first_fof_call;
  double time_first_fof_call;
  double delta_time_fof;

  /* Integer time of the next FoF black holes seeding call */
  integertime_t ti_next_fof;

  /* FOF information */
  int run_fof;
  int dump_catalogue_when_seeding;

  /* power spectrum information */
  double a_first_ps_output;
  double time_first_ps_output;
  double delta_time_ps;
  int ps_output_count;

  /* Output_List for the power spectrum */
  struct output_list *output_list_ps;

  /* Integer time of the next ps output */
  integertime_t ti_next_ps;

  /* Statistics information */
  double a_first_statistics;
  double time_first_statistics;
  double delta_time_statistics;

  /* Output_List for the stats */
  struct output_list *output_list_stats;

  /* Integer time of the next statistics dump */
  integertime_t ti_next_stats;

  /* File handle for the statistics */
  FILE *file_stats;

  /* File handle for the timesteps information */
  FILE *file_timesteps;

  /* File handle for the Radiative Transfer sub-cycling information */
  FILE *file_rt_subcycles;

  /* File handle for the SFH logger file */
  FILE *sfh_logger;

  /* The current step number. */
  int step;

  /* Data for the threads' barrier. */
  swift_barrier_t wait_barrier;
  swift_barrier_t run_barrier;

  /* ID of the node this engine lives on. */
  int nr_nodes, nodeID;

  /* Proxies for the other nodes in this simulation. */
  struct proxy *proxies;
  int nr_proxies, *proxy_ind;

  /* Tic/toc at the start/end of a step. */
  ticks tic_step, toc_step;

#ifdef WITH_MPI
  /* CPU times that the tasks used in the last step. */
  double usertime_last_step;
  double systime_last_step;

  /* Step of last repartition. */
  int last_repartition;

  /* Use synchronous redistributes. */
  int syncredist;

#endif

  /* Wallclock time of the last time-step */
  float wallclock_time;

  /* Are we in the process of restaring a simulation? */
  int restarting;

  /* Force the engine to rebuild? */
  int forcerebuild;

  /* Force the engine to repartition ? */
  int forcerepart;
  struct repartition *reparttype;

  /* The Continuous Simulation Data Stream (CSDS) */
  struct csds_writer *csds;

  /* How many steps have we done with the same set of tasks? */
  int tasks_age;

  /* Linked list for cell-task association. */
  struct link *links;
  size_t nr_links, size_links;

  /* Average number of tasks per cell. Used to estimate the sizes
   * of the various task arrays. Also the maximum from all ranks. */
  float tasks_per_cell;
  float tasks_per_cell_max;

  /* Average number of links per tasks. This number is used before
     the creation of communication tasks so needs to be large enough. */
  float links_per_tasks;

  /* Are we talkative ? */
  int verbose;

  /* Physical constants definition */
  const struct phys_const *physical_constants;

  /* The cosmological model */
  struct cosmology *cosmology;

  /* Properties of the hydro scheme */
  struct hydro_props *hydro_properties;

  /* Properties of the entropy floor */
  const struct entropy_floor_properties *entropy_floor;

  /* Properties of the star model */
  struct stars_props *stars_properties;

  /* Properties of the black hole model */
  const struct black_holes_props *black_holes_properties;

  /* Properties of the sink model */
  const struct sink_props *sink_properties;

  /* Properties of the neutrino model */
  const struct neutrino_props *neutrino_properties;

  /* The linear neutrino response */
  struct neutrino_response *neutrino_response;

  /* Properties of the self-gravity scheme */
  struct gravity_props *gravity_properties;

  /* The mesh used for long-range gravity forces */
  struct pm_mesh *mesh;

  /* Properties and pointers for the power spectrum */
  struct power_spectrum_data *power_data;

  /* Properties of external gravitational potential */
  const struct external_potential *external_potential;

  /* Properties of the hydrodynamics forcing terms */
  const struct forcing_terms *forcing_terms;

  /* Properties of the cooling scheme */
  struct cooling_function_data *cooling_func;

  /* Properties of the starformation law */
  const struct star_formation *star_formation;

  /* Properties of the sellar feedback model */
  struct feedback_props *feedback_props;

  /* Properties of the pressure floor scheme */
  struct pressure_floor_props *pressure_floor_props;

  /* Properties of the radiative transfer model */
  struct rt_props *rt_props;

  /* Properties of the chemistry model */
  const struct chemistry_global_data *chemistry;

  /* Properties used to compute the extra i/o fields */
  struct extra_io_properties *io_extra_props;

  /*! The FOF properties data. */
  struct fof_props *fof_properties;

  /* The (parsed) parameter file */
  struct swift_params *parameter_file;

  /* The output selection options */
  struct output_options *output_options;

  /* Temporary struct to hold a group of deferable properties (in MPI mode
   * these are reduced together, but may not be required just yet). */
  struct collectgroup1 collect_group1;

  /* Whether to dump restart files. */
  int restart_dump;

  /* Whether to save previous generation of restart files. */
  int restart_save;

  /* Whether to dump restart files after the last step. */
  int restart_onexit;

  /* Number of Lustre OSTs on the system to use as rank-based striping offset */
  int restart_lustre_OST_count;

  /* Do we free the foreign data before writing restart files? */
  int free_foreign_when_dumping_restart;

  /* Do we free the foreign data before rebuilding the tree? */
  int free_foreign_when_rebuilding;

  /* Name of the restart file directory. */
  const char *restart_dir;

  /* Name of the restart file. */
  const char *restart_file;

  /* Flag whether we should resubmit on this step? */
  int resubmit;

  /* Do we want to run the resubmission command once a run reaches the time
   * limit? */
  int resubmit_after_max_hours;

  /* What command should we run to resubmit at the end? */
  char resubmit_command[PARSER_MAX_LINE_SIZE];

  /* How often to check for the stop file and dump restarts and exit the
   * application. */
  int restart_stop_steps;

  /* Get the maximal wall-clock time of this run */
  float restart_max_hours_runtime;

  /* Ticks between restart dumps. */
  ticks restart_dt;

  /* Time after which next dump will occur. */
  ticks restart_next;

  /* Maximum number of tasks needed for restarting. */
  int restart_max_tasks;

  /* The globally agreed runtime, in hours. */
  float runtime;

  /* The locally accumulated deadtime. */
  double local_deadtime;

  /* The globally accumulated deadtime. */
  double global_deadtime;

  /* Time-integration mesh kick to apply to the particle velocities for
   * snapshots */
  float dt_kick_grav_mesh_for_io;

  /* Label of the run */
  char run_name[PARSER_MAX_LINE_SIZE];

  /* Has there been an stf this timestep? */
  char stf_this_timestep;

  /* Line of sight properties. */
  struct los_props *los_properties;

  /* Line of sight properties. */
  struct lightcone_array_props *lightcone_array_properties;

  /* Line of sight outputs information. */
  struct output_list *output_list_los;
  double a_first_los;
  double time_first_los;
  double delta_time_los;
  integertime_t ti_next_los;
  int los_output_count;

  /* Lightcone information */
  int flush_lightcone_maps;

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Run brute force checks only on steps when all gparts active? */
  int force_checks_only_all_active;

  /* Run brute force checks only during snapshot timesteps? */
  int force_checks_only_at_snapshots;

  /* Are all gparts active this timestep? */
  int all_gparts_active;

  /* Flag to tell brute force checks a snapshot was recently written. */
  int force_checks_snapshot_flag;
#endif
};

/* Function prototypes, engine.c. */
void engine_addlink(struct engine *e, struct link **l, struct task *t);
void engine_barrier(struct engine *e);
void engine_compute_next_snapshot_time(struct engine *e);
void engine_compute_next_stf_time(struct engine *e);
void engine_compute_next_fof_time(struct engine *e);
void engine_compute_next_statistics_time(struct engine *e);
void engine_compute_next_los_time(struct engine *e);
void engine_compute_next_ps_time(struct engine *e);
void engine_recompute_displacement_constraint(struct engine *e);
void engine_unskip(struct engine *e);
void engine_unskip_rt_sub_cycle(struct engine *e);
void engine_drift_all(struct engine *e, const int drift_mpoles);
void engine_drift_top_multipoles(struct engine *e);
void engine_reconstruct_multipoles(struct engine *e);
void engine_allocate_foreign_particles(struct engine *e, const int fof);
void engine_print_stats(struct engine *e);
void engine_io(struct engine *e);
void engine_io_check_snapshot_triggers(struct engine *e);
void engine_collect_end_of_step(struct engine *e, int apply);
void engine_collect_end_of_sub_cycle(struct engine *e);
void engine_dump_snapshot(struct engine *e, const int fof);
void engine_run_on_dump(struct engine *e);
void engine_init_output_lists(struct engine *e, struct swift_params *params,
                              const struct output_options *output_options);
void engine_init(
    struct engine *e, struct space *s, struct swift_params *params,
    struct output_options *output_options, long long Ngas, long long Ngparts,
    long long Nsinks, long long Nstars, long long Nblackholes,
    long long Nbackground_gparts, long long Nnuparts, int policy, int verbose,
    const struct unit_system *internal_units,
    const struct phys_const *physical_constants, struct cosmology *cosmo,
    struct hydro_props *hydro,
    const struct entropy_floor_properties *entropy_floor,
    struct gravity_props *gravity, struct stars_props *stars,
    const struct black_holes_props *black_holes, const struct sink_props *sinks,
    const struct neutrino_props *neutrinos,
    struct neutrino_response *neutrino_response,
    struct feedback_props *feedback,
    struct pressure_floor_props *pressure_floor, struct rt_props *rt,
    struct pm_mesh *mesh, struct power_spectrum_data *pow_data,
    const struct external_potential *potential,
    const struct forcing_terms *forcing_terms,
    struct cooling_function_data *cooling_func,
    const struct star_formation *starform,
    const struct chemistry_global_data *chemistry,
    struct extra_io_properties *io_extra_props,
    struct fof_props *fof_properties, struct los_props *los_properties,
    struct lightcone_array_props *lightcone_array_properties,
    struct ic_info *ics_metadata);
void engine_config(int restart, int fof, struct engine *e,
                   struct swift_params *params, int nr_nodes, int nodeID,
                   int nr_task_threads, int nr_pool_threads, int with_aff,
                   int verbose, const char *restart_dir,
                   const char *restart_file, struct repartition *reparttype);
void engine_launch(struct engine *e, const char *call);
int engine_prepare(struct engine *e);
void engine_run_rt_sub_cycles(struct engine *e);
void engine_init_particles(struct engine *e, int flag_entropy_ICs,
                           int clean_h_values);
int engine_step(struct engine *e);
void engine_split(struct engine *e, struct partition *initial_partition);
void engine_exchange_strays(struct engine *e, const size_t offset_parts,
                            const int *ind_part, size_t *Npart,
                            const size_t offset_gparts, const int *ind_gpart,
                            size_t *Ngpart, const size_t offset_sparts,
                            const int *ind_spart, size_t *Nspart,
                            const size_t offset_bparts, const int *ind_bpart,
                            size_t *Nbpart);
void engine_rebuild(struct engine *e, int redistributed, int clean_h_values);
void engine_repartition(struct engine *e);
void engine_repartition_trigger(struct engine *e);
void engine_makeproxies(struct engine *e);
void engine_redistribute(struct engine *e);
void engine_print_policy(struct engine *e);
int engine_is_done(struct engine *e);
void engine_pin(void);
void engine_unpin(void);
void engine_clean(struct engine *e, const int fof, const int restart);
int engine_estimate_nr_tasks(const struct engine *e);
void engine_print_task_counts(const struct engine *e);
void engine_fof(struct engine *e, const int dump_results,
                const int dump_debug_results, const int seed_black_holes,
                const int foreign_buffers_allocated);
void engine_activate_gpart_comms(struct engine *e);
void engine_activate_fof_attach_tasks(struct engine *e);

/* Function prototypes, engine_maketasks.c. */
void engine_maketasks(struct engine *e);

/* Function prototypes, engine_maketasks.c. */
void engine_make_fof_tasks(struct engine *e);

/* Function prototypes, engine_marktasks.c. */
int engine_marktasks(struct engine *e);

/* Function prototypes, engine_split_particles.c. */
void engine_split_gas_particles(struct engine *e);

#ifdef HAVE_SETAFFINITY
cpu_set_t *engine_entry_affinity(void);
#endif
void engine_numa_policies(int rank, int verbose);

/* Struct dump/restore support. */
void engine_struct_dump(struct engine *e, FILE *stream);
void engine_struct_restore(struct engine *e, FILE *stream);
int engine_dump_restarts(struct engine *e, int drifted_all, int force);

/* dev/debug */
void engine_dump_diagnostic_data(struct engine *e);

#endif /* SWIFT_ENGINE_H */
