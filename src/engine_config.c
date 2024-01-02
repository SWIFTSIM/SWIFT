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

/* System includes. */
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef HAVE_LIBNUMA
#include <numa.h>
#endif

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "fof.h"
#include "mpiuse.h"
#include "part.h"
#include "pressure_floor.h"
#include "proxy.h"
#include "rt.h"
#include "star_formation.h"
#include "star_formation_logger.h"
#include "stars_io.h"
#include "statistics.h"
#include "version.h"

extern int engine_max_parts_per_ghost;
extern int engine_max_sparts_per_ghost;
extern int engine_max_parts_per_cooling;

/**
 * @brief dump diagnostic data on tasks, memuse, mpiuse, queues.
 *
 * @param e the #engine
 */
void engine_dump_diagnostic_data(struct engine *e) {
  /* OK, do our work. */
  message("Dumping engine tasks in step: %d", e->step);
  task_dump_active(e);

#ifdef SWIFT_MEMUSE_REPORTS
  /* Dump the currently logged memory. */
  message("Dumping memory use report");
  memuse_log_dump_error(e->nodeID);
#endif

#if defined(SWIFT_MPIUSE_REPORTS) && defined(WITH_MPI)
  /* Dump the MPI interactions in the step. */
  mpiuse_log_dump_error(e->nodeID);
#endif

  /* Add more interesting diagnostics. */
  scheduler_dump_queues(e);
}

/* Particle cache size. */
#define CACHE_SIZE 512

#ifdef SWIFT_DUMPER_THREAD
/**
 * @brief dumper thread action, checks got the existence of the .dump file
 * every 5 seconds and does the dump if found.
 *
 * @param p the #engine
 */
static void *engine_dumper_poll(void *p) {
  struct engine *e = (struct engine *)p;

#ifdef WITH_MPI
  char dumpfile[10];
  snprintf(dumpfile, sizeof(dumpfile), ".dump.%d", e->nodeID);
#else
  const char *dumpfile = ".dump";
#endif

  while (1) {
    if (access(dumpfile, F_OK) == 0) {

      engine_dump_diagnostic_data(e);

      /* Delete the file. */
      unlink(dumpfile);
      message("Dumping completed");
      fflush(stdout);
    }

    /* Take a breath. */
    sleep(5);
  }
  return NULL;
}
#endif /* SWIFT_DUMPER_THREAD */

#ifdef SWIFT_DUMPER_THREAD
/**
 * @brief creates the dumper thread.
 *
 * This watches for the creation of a ".dump" file in the current directory
 * and if found dumps the current state of the tasks and memory use (if also
 * configured).
 *
 * @param e the #engine
 *
 */
static void engine_dumper_init(struct engine *e) {
  pthread_t dumper;

#ifdef WITH_MPI
  char dumpfile[10];
  snprintf(dumpfile, sizeof(dumpfile), ".dump.%d", e->nodeID);
#else
  const char *dumpfile = ".dump";
#endif

  /* Make sure the .dump file is not present, that is bad when starting up. */
  struct stat buf;
  if (stat(dumpfile, &buf) == 0) unlink(dumpfile);

  /* Thread does not exit, so nothing to do but create it. */
  pthread_create(&dumper, NULL, &engine_dumper_poll, e);
}
#endif /* SWIFT_DUMPER_THREAD */

/**
 * @brief configure an engine with the given number of threads, queues
 *        and core affinity. Also initialises the scheduler and opens various
 *        output files, computes the next timestep and initialises the
 *        threadpool.
 *
 * Assumes the engine is correctly initialised i.e. is restored from a restart
 * file or has been setup by engine_init(). When restarting any output log
 * files are positioned so that further output is appended. Note that
 * parameters are not read from the engine, just the parameter file, this
 * allows values derived in this function to be changed between runs.
 * When not restarting params should be the same as given to engine_init().
 *
 * @param restart true when restarting the application.
 * @param fof true when starting a stand-alone FOF call.
 * @param e The #engine.
 * @param params The parsed parameter file.
 * @param nr_nodes The number of MPI ranks.
 * @param nodeID The MPI rank of this node.
 * @param nr_task_threads The number of engine threads per MPI rank.
 * @param nr_pool_threads The number of threadpool threads per MPI rank.
 * @param with_aff use processor affinity, if supported.
 * @param verbose Is this #engine talkative ?
 * @param restart_file The name of our restart file.
 * @param reparttype What type of repartition algorithm are we using.
 */
void engine_config(int restart, int fof, struct engine *e,
                   struct swift_params *params, int nr_nodes, int nodeID,
                   int nr_task_threads, int nr_pool_threads, int with_aff,
                   int verbose, const char *restart_dir,
                   const char *restart_file, struct repartition *reparttype) {

  struct clocks_time tic, toc;
  if (nodeID == 0) clocks_gettime(&tic);

  /* Store the values and initialise global fields. */
  e->nodeID = nodeID;
  e->nr_threads = nr_task_threads;
  e->nr_pool_threads = nr_pool_threads;
  e->nr_nodes = nr_nodes;
  e->proxy_ind = NULL;
  e->nr_proxies = 0;
  e->forcerebuild = 1;
  e->forcerepart = 0;
  e->restarting = restart;
  e->step_props = engine_step_prop_none;
  e->links = NULL;
  e->nr_links = 0;
  e->file_stats = NULL;
  e->file_timesteps = NULL;
  e->file_rt_subcycles = NULL;
  e->sfh_logger = NULL;
  e->verbose = verbose;
  e->wallclock_time = 0.f;
  e->restart_dump = 0;
  e->restart_dir = restart_dir;
  e->restart_file = restart_file;
  e->resubmit = 0;
  e->restart_next = 0;
  e->restart_dt = 0;
  e->run_fof = 0;

  /* Seed rand(). */
  srand(clocks_random_seed());

  /* Allow repartitioning to be changed between restarts. On restart this is
   * already allocated and freed on exit, so we need to copy over. */
#ifdef WITH_MPI
  if (restart) {
    int *celllist = e->reparttype->celllist;
    int ncelllist = e->reparttype->ncelllist;
    memcpy(e->reparttype, reparttype, sizeof(struct repartition));
    e->reparttype->celllist = celllist;
    e->reparttype->ncelllist = ncelllist;
  } else {
    e->reparttype = reparttype;
  }
#endif

  if (restart && fof) {
    error(
        "Can't configure the engine to be a stand-alone FOF and restarting "
        "from a check-point at the same time!");
  }

  /* Welcome message */
  if (e->nodeID == 0) message("Running simulation '%s'.", e->run_name);

  /* Check-pointing properties */

  e->restart_stop_steps =
      parser_get_opt_param_int(params, "Restarts:stop_steps", 100);

  e->restart_max_hours_runtime =
      parser_get_opt_param_float(params, "Restarts:max_run_time", FLT_MAX);

  e->resubmit_after_max_hours =
      parser_get_opt_param_int(params, "Restarts:resubmit_on_exit", 0);

  if (e->resubmit_after_max_hours)
    parser_get_param_string(params, "Restarts:resubmit_command",
                            e->resubmit_command);

  /* Get the number of queues */
  int nr_queues =
      parser_get_opt_param_int(params, "Scheduler:nr_queues", e->nr_threads);
  if (nr_queues <= 0) nr_queues = e->nr_threads;
  if (nr_queues != nr_task_threads)
    message("Number of task queues set to %d", nr_queues);
  e->s->nr_queues = nr_queues;

  /* Get the frequency of the dependency graph dumping */
  e->sched.frequency_dependency = parser_get_opt_param_int(
      params, "Scheduler:dependency_graph_frequency", 0);
  if (e->sched.frequency_dependency < 0) {
    error("Scheduler:dependency_graph_frequency should be >= 0");
  }
  /* Get cellID for extra dependency graph dumps of specific cell */
  e->sched.dependency_graph_cellID = parser_get_opt_param_longlong(
      params, "Scheduler:dependency_graph_cell", 0LL);

  /* Get the frequency of the task level dumping */
  e->sched.frequency_task_levels = parser_get_opt_param_int(
      params, "Scheduler:task_level_output_frequency", 0);
  if (e->sched.frequency_task_levels < 0) {
    error("Scheduler:task_level_output_frequency should be >= 0");
  }

#if defined(SWIFT_DEBUG_CHECKS)
  e->sched.deadlock_waiting_time_ms = parser_get_opt_param_float(
      params, "Scheduler:deadlock_waiting_time_s", -1.f);
  /* User provides parameter in s. We want it in ms. */
  e->sched.deadlock_waiting_time_ms *= 1000.f;
#endif

/* Deal with affinity. For now, just figure out the number of cores. */
#if defined(HAVE_SETAFFINITY)
  const int nr_cores = sysconf(_SC_NPROCESSORS_ONLN);
  cpu_set_t *entry_affinity = engine_entry_affinity();
  const int nr_affinity_cores = CPU_COUNT(entry_affinity);

  if (nr_cores > CPU_SETSIZE) /* Unlikely, except on e.g. SGI UV. */
    error("must allocate dynamic cpu_set_t (too many cores per node)");

  if (verbose && with_aff) {
    char *buf = (char *)malloc((nr_cores + 1) * sizeof(char));
    buf[nr_cores] = '\0';
    for (int j = 0; j < nr_cores; ++j) {
      /* Reversed bit order from convention, but same as e.g. Intel MPI's
       * I_MPI_PIN_DOMAIN explicit mask: left-to-right, LSB-to-MSB. */
      buf[j] = CPU_ISSET(j, entry_affinity) ? '1' : '0';
    }
    message("Affinity at entry: %s", buf);
    free(buf);
  }

  int *cpuid = NULL;
  cpu_set_t cpuset;

  if (with_aff) {

    cpuid = (int *)malloc(nr_affinity_cores * sizeof(int));

    int skip = 0;
    for (int k = 0; k < nr_affinity_cores; k++) {
      int c;
      for (c = skip; c < CPU_SETSIZE && !CPU_ISSET(c, entry_affinity); ++c)
        ;
      cpuid[k] = c;
      skip = c + 1;
    }

#if defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
    if ((e->policy & engine_policy_cputight) != engine_policy_cputight) {

      if (numa_available() >= 0) {
        if (nodeID == 0) message("prefer NUMA-distant CPUs");

        /* Get list of numa nodes of all available cores. */
        int *nodes = (int *)malloc(nr_affinity_cores * sizeof(int));
        int nnodes = 0;
        for (int i = 0; i < nr_affinity_cores; i++) {
          nodes[i] = numa_node_of_cpu(cpuid[i]);
          if (nodes[i] > nnodes) nnodes = nodes[i];
        }
        nnodes += 1;

        /* Count cores per node. */
        int *core_counts = (int *)malloc(nnodes * sizeof(int));
        for (int i = 0; i < nr_affinity_cores; i++) {
          core_counts[nodes[i]] = 0;
        }
        for (int i = 0; i < nr_affinity_cores; i++) {
          core_counts[nodes[i]] += 1;
        }

        /* Index cores within each node. */
        int *core_indices = (int *)malloc(nr_affinity_cores * sizeof(int));
        for (int i = nr_affinity_cores - 1; i >= 0; i--) {
          core_indices[i] = core_counts[nodes[i]];
          core_counts[nodes[i]] -= 1;
        }

        /* Now sort so that we pick adjacent cpuids from different nodes
         * by sorting internal node core indices. */
        int done = 0;
        while (!done) {
          done = 1;
          for (int i = 1; i < nr_affinity_cores; i++) {
            if (core_indices[i] < core_indices[i - 1]) {
              int t = cpuid[i - 1];
              cpuid[i - 1] = cpuid[i];
              cpuid[i] = t;

              t = core_indices[i - 1];
              core_indices[i - 1] = core_indices[i];
              core_indices[i] = t;
              done = 0;
            }
          }
        }

        free(nodes);
        free(core_counts);
        free(core_indices);
      }
    }
#endif
  } else {
    if (nodeID == 0) message("no processor affinity used");

  } /* with_aff */

  /* Avoid (unexpected) interference between engine and runner threads. We can
   * do this once we've made at least one call to engine_entry_affinity and
   * maybe numa_node_of_cpu(sched_getcpu()), even if the engine isn't already
   * pinned. */
  if (with_aff) engine_unpin();
#endif

  if (with_aff && nodeID == 0) {
#ifdef HAVE_SETAFFINITY
#ifdef WITH_MPI
    printf("[%04i] %s engine_init: cpu map is [ ", nodeID,
           clocks_get_timesincestart());
#else
    printf("%s engine_init: cpu map is [ ", clocks_get_timesincestart());
#endif
    for (int i = 0; i < nr_affinity_cores; i++) printf("%i ", cpuid[i]);
    printf("].\n");
#endif
  }

  /* Are we doing stuff in parallel? */
  if (nr_nodes > 1) {
#ifndef WITH_MPI
    error("SWIFT was not compiled with MPI support.");
#else

    /* Make sure the corresponding policy is set and make space for the proxies
     */
    e->policy |= engine_policy_mpi;
    if ((e->proxies = (struct proxy *)calloc(sizeof(struct proxy),
                                             engine_maxproxies)) == NULL)
      error("Failed to allocate memory for proxies.");
    e->nr_proxies = 0;

    /* Use synchronous MPI sends and receives when redistributing. */
    e->syncredist =
        parser_get_opt_param_int(params, "DomainDecomposition:synchronous", 0);

    /* Collect the hostname of each rank into a file */

    const int hostname_buffer_length = 256;
    char my_hostname[256] = {0};
    sprintf(my_hostname, "%s", hostname());

    char *hostnames = NULL;

    if (nodeID == 0) {
      hostnames =
          (char *)calloc(nr_nodes * hostname_buffer_length, sizeof(char));
      if (hostnames == NULL)
        error("Failed to allocate memory for hostname list");
    }

    MPI_Gather(my_hostname, hostname_buffer_length, MPI_BYTE, hostnames,
               hostname_buffer_length, MPI_BYTE, 0, MPI_COMM_WORLD);

    if (nodeID == 0) {
      FILE *ranklog = NULL;
      if (restart) {
        ranklog = fopen("rank_hostname.log", "a");
        if (ranklog == NULL)
          error("Could not create file 'rank_hostname.log'.");
      } else {
        ranklog = fopen("rank_hostname.log", "w");
        if (ranklog == NULL)
          error("Could not open file 'rank_hostname.log' for writing.");
      }

      /* Write the header every restart-cycle. It does not hurt. */
      fprintf(ranklog, "# step rank hostname\n");

      for (int i = 0; i < nr_nodes; ++i)
        fprintf(ranklog, "%d %d %s\n", e->step, i,
                hostnames + hostname_buffer_length * i);

      fclose(ranklog);
      free(hostnames);
    }
#endif
  }

  /* Open some global files */
  if (!fof && e->nodeID == 0) {

    /* When restarting append to these files. */
    const char *mode;
    if (restart)
      mode = "a";
    else
      mode = "w";

    char energyfileName[200] = "";
    parser_get_opt_param_string(params, "Statistics:energy_file_name",
                                energyfileName,
                                engine_default_energy_file_name);
    sprintf(energyfileName + strlen(energyfileName), ".txt");
    e->file_stats = fopen(energyfileName, mode);
    if (e->file_stats == NULL)
      error("Could not open the file '%s' with mode '%s'.", energyfileName,
            mode);

    if (!restart)
      stats_write_file_header(e->file_stats, e->internal_units,
                              e->physical_constants);

    char timestepsfileName[200] = "";
    parser_get_opt_param_string(params, "Statistics:timestep_file_name",
                                timestepsfileName,
                                engine_default_timesteps_file_name);

    sprintf(timestepsfileName + strlen(timestepsfileName), ".txt");
    e->file_timesteps = fopen(timestepsfileName, mode);
    if (e->file_timesteps == NULL)
      error("Could not open the file '%s' with mode '%s'.", timestepsfileName,
            mode);

#ifndef RT_NONE
    char rtSubcyclesFileName[200] = "";
    parser_get_opt_param_string(params, "Statistics:rt_subcycles_file_name",
                                rtSubcyclesFileName,
                                engine_default_rt_subcycles_file_name);
    sprintf(rtSubcyclesFileName + strlen(rtSubcyclesFileName), ".txt");
    e->file_rt_subcycles = fopen(rtSubcyclesFileName, mode);
    if (e->file_rt_subcycles == NULL)
      error("Could not open the file '%s' with mode '%s'.", rtSubcyclesFileName,
            mode);
#endif

    if (!restart) {
      fprintf(
          e->file_timesteps,
          "# Host: %s\n# Branch: %s\n# Revision: %s\n# Compiler: %s, "
          "Version: %s \n# "
          "Number of threads: %d\n# Number of MPI ranks: %d\n# Hydrodynamic "
          "scheme: %s\n# Hydrodynamic kernel: %s\n# No. of neighbours: %.2f "
          "+/- %.4f\n# Eta: %f\n# Config: %s\n# CFLAGS: %s\n",
          hostname(), git_branch(), git_revision(), compiler_name(),
          compiler_version(), e->nr_threads, e->nr_nodes, SPH_IMPLEMENTATION,
          kernel_name, e->hydro_properties->target_neighbours,
          e->hydro_properties->delta_neighbours,
          e->hydro_properties->eta_neighbours, configuration_options(),
          compilation_cflags());

      fprintf(
          e->file_timesteps,
          "# Step Properties: Rebuild=%d, Redistribute=%d, Repartition=%d, "
          "Statistics=%d, Snapshot=%d, Restarts=%d STF=%d, FOF=%d, mesh=%d\n",
          engine_step_prop_rebuild, engine_step_prop_redistribute,
          engine_step_prop_repartition, engine_step_prop_statistics,
          engine_step_prop_snapshot, engine_step_prop_restarts,
          engine_step_prop_stf, engine_step_prop_fof, engine_step_prop_mesh);

      fprintf(
          e->file_timesteps,
          "# %6s %14s %12s %12s %14s %9s %12s %12s %12s %12s %12s %16s "
          "[%s] %6s %12s [%s] "
          "%16s [%2s] %16s [%2s] %16s [%2s] %16s [%2s] %16s [%2s] %16s [%2s] "
          "%16s [%2s] %16s [%2s] %16s [%2s] %16s [%2s] %16s [%2s] %16s [%2s]\n",
          "Step", "Time", "Scale-factor", "Redshift", "Time-step", "Time-bins",
          "Updates", "g-Updates", "s-Updates", "Sink-Updates", "b-Updates",
          "Wall-clock time", clocks_getunit(), "Props", "Dead time",
          clocks_getunit(), "Drifts", clocks_getunit(), "Sorts",
          clocks_getunit(), "Hyrdo", clocks_getunit(), "Hydro Ghost",
          clocks_getunit(), "Hydro Density", clocks_getunit(), "Gravity",
          clocks_getunit(), "Feedback", clocks_getunit(), "Time Integration",
          clocks_getunit(), "MPI", clocks_getunit(), "RT", clocks_getunit(),
          "Tchem", clocks_getunit(), "Others", clocks_getunit());
      fflush(e->file_timesteps);

#ifndef RT_NONE
      fprintf(
          e->file_rt_subcycles,
          "# Host: %s\n# Branch: %s\n# Revision: %s\n# Compiler: %s, "
          "Version: %s \n# "
          "Number of threads: %d\n# Number of MPI ranks: %d\n# Hydrodynamic "
          "scheme: %s\n# Hydrodynamic kernel: %s\n# No. of neighbours: %.2f "
          "+/- %.4f\n# Eta: %f\n# Radiative Transfer Scheme: %s\n# Max Number "
          "RT sub-cycles: %d\n# Config: %s\n# CFLAGS: %s\n",
          hostname(), git_branch(), git_revision(), compiler_name(),
          compiler_version(), e->nr_threads, e->nr_nodes, SPH_IMPLEMENTATION,
          kernel_name, e->hydro_properties->target_neighbours,
          e->hydro_properties->delta_neighbours,
          e->hydro_properties->eta_neighbours, RT_IMPLEMENTATION,
          e->max_nr_rt_subcycles, configuration_options(),
          compilation_cflags());

      fprintf(
          e->file_rt_subcycles,
          "# Step Properties: Rebuild=%d, Redistribute=%d, Repartition=%d, "
          "Statistics=%d, Snapshot=%d, Restarts=%d STF=%d, FOF=%d, mesh=%d\n",
          engine_step_prop_rebuild, engine_step_prop_redistribute,
          engine_step_prop_repartition, engine_step_prop_statistics,
          engine_step_prop_snapshot, engine_step_prop_restarts,
          engine_step_prop_stf, engine_step_prop_fof, engine_step_prop_mesh);

      fprintf(e->file_rt_subcycles,
              "# Note: Sub-cycle=0 is performed during the regular SWIFT step, "
              "alongside hydro, gravity etc.\n");
      fprintf(e->file_rt_subcycles,
              "#       For this reason, the wall-clock time and dead time is "
              "not available for it, and is written as -1.\n");

      fprintf(e->file_rt_subcycles,
              "# %6s %9s %14s %12s %12s %14s %9s %12s %16s [%s] %12s [%s]\n",
              "Step", "Sub-cycle", "Time", "Scale-factor", "Redshift",
              "Time-step", "Time-bins", "RT-Updates", "Wall-clock time",
              clocks_getunit(), "Dead time", clocks_getunit());
      fflush(e->file_rt_subcycles);
#endif  // compiled with RT
    }

    /* Initialize the SFH logger if running with star formation */
    if (e->policy & engine_policy_star_formation) {
      e->sfh_logger = fopen("SFR.txt", mode);
      if (e->sfh_logger == NULL)
        error("Could not open the file 'SFR.txt' with mode '%s'.", mode);

      if (!restart) {
        star_formation_logger_init_log_file(e->sfh_logger, e->internal_units,
                                            e->physical_constants);
        fflush(e->sfh_logger);
      }
    }
  }

  /* Print policy */
  engine_print_policy(e);

  if (!fof) {

    /* Print information about the hydro scheme */
    if (e->policy & engine_policy_hydro) {
      if (e->nodeID == 0) hydro_props_print(e->hydro_properties);
      if (e->nodeID == 0) pressure_floor_print(e->pressure_floor_props);
      if (e->nodeID == 0) entropy_floor_print(e->entropy_floor);
    }

    /* Print information about the gravity scheme */
    if (e->policy & engine_policy_self_gravity)
      if (e->nodeID == 0) gravity_props_print(e->gravity_properties);

    /* Print information about the stellar scheme */
    if (e->policy & engine_policy_stars)
      if (e->nodeID == 0) stars_props_print(e->stars_properties);

    /* Print information about the RT scheme */
    if (e->policy & engine_policy_rt) {
      rt_props_print(e->rt_props);
      if (e->nodeID == 0) {
        if (e->max_nr_rt_subcycles <= 1)
          message("WARNING: running without RT sub-cycling.");
        else {
          /* Make sure max_nr_rt_subcycles is an acceptable power of 2 */
          timebin_t power_subcycles = 0;
          while ((e->max_nr_rt_subcycles > (1 << power_subcycles)) &&
                 power_subcycles < num_time_bins)
            ++power_subcycles;
          if (power_subcycles == num_time_bins)
            error("TimeIntegration:max_nr_rt_subcycles=%d too big",
                  e->max_nr_rt_subcycles);
          if ((1 << power_subcycles) > e->max_nr_rt_subcycles)
            error("TimeIntegration:max_nr_rt_subcycles=%d not a power of 2",
                  e->max_nr_rt_subcycles);
          message("Running up to %d RT sub-cycles per hydro step.",
                  e->max_nr_rt_subcycles);
        }
      }
    }

    /* Check we have sensible time bounds */
    if (e->time_begin >= e->time_end)
      error(
          "Final simulation time (t_end = %e) must be larger than the start "
          "time (t_beg = %e)",
          e->time_end, e->time_begin);

    /* Check we have sensible time-step values */
    if (e->dt_min > e->dt_max)
      error(
          "Minimal time-step size (%e) must be smaller than maximal time-step "
          "size (%e)",
          e->dt_min, e->dt_max);

    /* Info about time-steps */
    if (e->nodeID == 0) {
      message("Absolute minimal timestep size: %e", e->time_base);

      float dt_min = e->time_end - e->time_begin;
      while (dt_min > e->dt_min) dt_min /= 2.f;

      message("Minimal timestep size (on time-line): %e", dt_min);

      float dt_max = e->time_end - e->time_begin;
      while (dt_max > e->dt_max) dt_max /= 2.f;

      message("Maximal timestep size (on time-line): %e", dt_max);
    }

    if (e->dt_min < e->time_base && e->nodeID == 0)
      error(
          "Minimal time-step size smaller than the absolute possible minimum "
          "dt=%e",
          e->time_base);

    if (!(e->policy & engine_policy_cosmology))
      if (e->dt_max > (e->time_end - e->time_begin) && e->nodeID == 0)
        error("Maximal time-step size larger than the simulation run time t=%e",
              e->time_end - e->time_begin);

    /* Read (or re-read the list of outputs */
    engine_init_output_lists(e, params, e->output_options);

    /* Check whether output quantities make sense */
    if (e->policy & engine_policy_cosmology) {

      if (e->delta_time_snapshot <= 1.)
        error("Time between snapshots (%e) must be > 1.",
              e->delta_time_snapshot);

      if (e->delta_time_statistics <= 1.)
        error("Time between statistics (%e) must be > 1.",
              e->delta_time_statistics);

      if (e->a_first_snapshot < e->cosmology->a_begin)
        error(
            "Scale-factor of first snapshot (%e) must be after the simulation "
            "start a=%e.",
            e->a_first_snapshot, e->cosmology->a_begin);

      if (e->a_first_statistics < e->cosmology->a_begin)
        error(
            "Scale-factor of first stats output (%e) must be after the "
            "simulation start a=%e.",
            e->a_first_statistics, e->cosmology->a_begin);

      if (e->policy & engine_policy_structure_finding) {

        if (e->delta_time_stf == -1. && !e->snapshot_invoke_stf)
          error("A value for `StructureFinding:delta_time` must be specified");

        if (e->a_first_stf_output < e->cosmology->a_begin)
          error(
              "Scale-factor of first stf output (%e) must be after the "
              "simulation start a=%e.",
              e->a_first_stf_output, e->cosmology->a_begin);
      }

      if (e->policy & engine_policy_fof &&
          e->fof_properties->seed_black_holes_enabled) {

        if (e->delta_time_fof <= 1.)
          error("Time between FOF (%e) must be > 1.", e->delta_time_fof);

        if (e->a_first_fof_call < e->cosmology->a_begin)
          error(
              "Scale-factor of first fof call (%e) must be after the "
              "simulation start a=%e.",
              e->a_first_fof_call, e->cosmology->a_begin);
      }

    } else {

      if (e->delta_time_snapshot <= 0.)
        error("Time between snapshots (%e) must be positive.",
              e->delta_time_snapshot);

      if (e->delta_time_statistics <= 0.)
        error("Time between statistics (%e) must be positive.",
              e->delta_time_statistics);

      /* Find the time of the first output */
      if (e->time_first_snapshot < e->time_begin)
        error(
            "Time of first snapshot (%e) must be after the simulation start "
            "t=%e.",
            e->time_first_snapshot, e->time_begin);

      if (e->time_first_statistics < e->time_begin)
        error(
            "Time of first stats output (%e) must be after the simulation "
            "start t=%e.",
            e->time_first_statistics, e->time_begin);

      if (e->policy & engine_policy_structure_finding) {

        if (e->delta_time_stf == -1. && !e->snapshot_invoke_stf)
          error("A value for `StructureFinding:delta_time` must be specified");

        if (e->delta_time_stf <= 0. && e->delta_time_stf != -1.)
          error("Time between STF (%e) must be positive.", e->delta_time_stf);

        if (e->time_first_stf_output < e->time_begin)
          error(
              "Time of first STF (%e) must be after the simulation start t=%e.",
              e->time_first_stf_output, e->time_begin);
      }
    }

    /* Get the total mass */
    e->total_mass = 0.;
    for (size_t i = 0; i < e->s->nr_gparts; ++i)
      e->total_mass += e->s->gparts[i].mass;

/* Reduce the total mass */
#ifdef WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, &e->total_mass, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

#if defined(WITH_CSDS)
    if ((e->policy & engine_policy_csds) && e->nodeID == 0)
      message(
          "WARNING: There is currently no way of predicting the output "
          "size, please use the CSDS carefully");
#endif

    /* Find the time of the first snapshot output */
    engine_compute_next_snapshot_time(e);

    /* Find the time of the first statistics output */
    engine_compute_next_statistics_time(e);

    /* Find the time of the first line of sight output */
    if (e->policy & engine_policy_line_of_sight) {
      engine_compute_next_los_time(e);
    }

    /* Find the time of the first stf output */
    if (e->policy & engine_policy_structure_finding) {
      engine_compute_next_stf_time(e);
    }

    /* Find the time of the first stf output */
    if (e->policy & engine_policy_fof &&
        e->fof_properties->seed_black_holes_enabled) {
      engine_compute_next_fof_time(e);
    }

    /* Check that we are invoking VELOCIraptor only if we have it */
    if (e->snapshot_invoke_stf &&
        !(e->policy & engine_policy_structure_finding)) {
      error(
          "Invoking VELOCIraptor after snapshots but structure finding wasn't "
          "activated at runtime (Use --velociraptor).");
    }

    /* Whether restarts are enabled. Yes by default. Can be changed on restart.
     */
    e->restart_dump = parser_get_opt_param_int(params, "Restarts:enable", 1);

    /* Whether to save backup copies of the previous restart files. */
    e->restart_save = parser_get_opt_param_int(params, "Restarts:save", 1);

    /* Whether restarts should be dumped on exit. Not by default. Can be changed
     * on restart. */
    e->restart_onexit = parser_get_opt_param_int(params, "Restarts:onexit", 0);

    /* Read the number of Lustre OSTs to distribute the restart files over */
    e->restart_lustre_OST_count =
        parser_get_opt_param_int(params, "Restarts:lustre_OST_count", 0);

    /* Hours between restart dumps. Can be changed on restart. */
    float dhours =
        parser_get_opt_param_float(params, "Restarts:delta_hours", 5.0f);
    if (e->nodeID == 0) {
      if (e->restart_dump)
        message("Restarts will be dumped every %f hours", dhours);
      else
        message("WARNING: restarts will not be dumped");

      if (e->verbose && e->restart_onexit)
        message("Restarts will be dumped after the final step");
    }

    /* Internally we use ticks, so convert into a delta ticks. Assumes we can
     * convert from ticks into milliseconds. */
    e->restart_dt = clocks_to_ticks(dhours * 60.0 * 60.0 * 1000.0);

    /* The first dump will happen no sooner than restart_dt ticks in the
     * future. */
    e->restart_next = getticks() + e->restart_dt;
  }

/* Construct types for MPI communications */
#ifdef WITH_MPI
  part_create_mpi_types();
  multipole_create_mpi_types();
  stats_create_mpi_type();
  proxy_create_mpi_type();
  task_create_mpi_comms();
#ifdef WITH_FOF
  fof_create_mpi_types();
#endif /* WITH_FOF */
#endif /* WITH_MPI */

  if (!fof) {

    /* Initialise the collection group. */
    collectgroup_init();
  }

  /* Initialize the threadpool. */
  threadpool_init(&e->threadpool, nr_pool_threads);
  if (e->nodeID == 0)
    message("Using %d threads in the thread-pool", nr_pool_threads);

  /* Cells per thread buffer. */
  e->s->cells_sub =
      (struct cell **)calloc(nr_pool_threads + 1, sizeof(struct cell *));
  e->s->multipoles_sub = (struct gravity_tensors **)calloc(
      nr_pool_threads + 1, sizeof(struct gravity_tensors *));

  /* First of all, init the barrier and lock it. */
  if (swift_barrier_init(&e->wait_barrier, NULL, e->nr_threads + 1) != 0 ||
      swift_barrier_init(&e->run_barrier, NULL, e->nr_threads + 1) != 0)
    error("Failed to initialize barrier.");

  /* Expected average for tasks per cell. If set to zero we use a heuristic
   * guess based on the numbers of cells and how many tasks per cell we expect.
   * On restart this number cannot be estimated (no cells yet), so we recover
   * from the end of the dumped run. Can be changed on restart. */
  e->tasks_per_cell =
      parser_get_opt_param_float(params, "Scheduler:tasks_per_cell", 0.0);
  e->tasks_per_cell_max = 0.0f;

  float maxtasks = 0;
  if (restart)
    maxtasks = e->restart_max_tasks;
  else
    maxtasks = engine_estimate_nr_tasks(e);

  /* Estimated number of links per tasks */
  e->links_per_tasks =
      parser_get_opt_param_float(params, "Scheduler:links_per_tasks", 25.);

  /* Init the scheduler. */
  scheduler_init(&e->sched, e->s, maxtasks, nr_queues,
                 (e->policy & scheduler_flag_steal), e->nodeID, &e->threadpool);

  /* Maximum size of MPI task messages, in KB, that should not be buffered,
   * that is sent using MPI_Issend, not MPI_Isend. 4Mb by default. Can be
   * changed on restart.
   */
  e->sched.mpi_message_limit =
      parser_get_opt_param_int(params, "Scheduler:mpi_message_limit", 4) * 1024;

  if (restart) {

    /* Overwrite the constants for the scheduler */
    space_maxsize = parser_get_opt_param_int(params, "Scheduler:cell_max_size",
                                             space_maxsize);
    space_subsize_pair_hydro = parser_get_opt_param_int(
        params, "Scheduler:cell_sub_size_pair_hydro", space_subsize_pair_hydro);
    space_subsize_self_hydro = parser_get_opt_param_int(
        params, "Scheduler:cell_sub_size_self_hydro", space_subsize_self_hydro);
    space_subsize_pair_stars = parser_get_opt_param_int(
        params, "Scheduler:cell_sub_size_pair_stars", space_subsize_pair_stars);
    space_subsize_self_stars = parser_get_opt_param_int(
        params, "Scheduler:cell_sub_size_self_stars", space_subsize_self_stars);
    space_subsize_pair_grav = parser_get_opt_param_int(
        params, "Scheduler:cell_sub_size_pair_grav", space_subsize_pair_grav);
    space_subsize_self_grav = parser_get_opt_param_int(
        params, "Scheduler:cell_sub_size_self_grav", space_subsize_self_grav);
    space_splitsize = parser_get_opt_param_int(
        params, "Scheduler:cell_split_size", space_splitsize);
    space_subdepth_diff_grav =
        parser_get_opt_param_int(params, "Scheduler:cell_subdepth_diff_grav",
                                 space_subdepth_diff_grav_default);
    space_extra_parts = parser_get_opt_param_int(
        params, "Scheduler:cell_extra_parts", space_extra_parts_default);
    space_extra_sparts = parser_get_opt_param_int(
        params, "Scheduler:cell_extra_sparts", space_extra_sparts_default);
    space_extra_gparts = parser_get_opt_param_int(
        params, "Scheduler:cell_extra_gparts", space_extra_gparts_default);
    space_extra_bparts = parser_get_opt_param_int(
        params, "Scheduler:cell_extra_bparts", space_extra_bparts_default);

    /* Do we want any spare particles for on the fly creation?
       This condition should be the same than in space.c */
    if (!(e->policy & engine_policy_star_formation ||
          e->policy & engine_policy_sinks) ||
        !swift_star_formation_model_creates_stars) {
      space_extra_sparts = 0;
      space_extra_gparts = 0;
      space_extra_sinks = 0;
    }

    engine_max_parts_per_ghost =
        parser_get_opt_param_int(params, "Scheduler:engine_max_parts_per_ghost",
                                 engine_max_parts_per_ghost);
    engine_max_sparts_per_ghost = parser_get_opt_param_int(
        params, "Scheduler:engine_max_sparts_per_ghost",
        engine_max_sparts_per_ghost);

    engine_max_parts_per_cooling = parser_get_opt_param_int(
        params, "Scheduler:engine_max_parts_per_cooling",
        engine_max_parts_per_cooling);
  }

  /* Allocate and init the threads. */
  if (swift_memalign("runners", (void **)&e->runners, SWIFT_CACHE_ALIGNMENT,
                     e->nr_threads * sizeof(struct runner)) != 0)
    error("Failed to allocate threads array.");

  for (int k = 0; k < e->nr_threads; k++) {
    e->runners[k].id = k;
    e->runners[k].e = e;
    if (pthread_create(&e->runners[k].thread, NULL, &runner_main,
                       &e->runners[k]) != 0)
      error("Failed to create runner thread.");

    /* Try to pin the runner to a given core */
    if (with_aff &&
        (e->policy & engine_policy_setaffinity) == engine_policy_setaffinity) {
#if defined(HAVE_SETAFFINITY)

      /* Set a reasonable queue ID. */
      int coreid = k % nr_affinity_cores;
      e->runners[k].cpuid = cpuid[coreid];

      if (nr_queues < e->nr_threads)
        e->runners[k].qid = cpuid[coreid] * nr_queues / nr_affinity_cores;
      else
        e->runners[k].qid = k;

      /* Set the cpu mask to zero | e->id. */
      CPU_ZERO(&cpuset);
      CPU_SET(cpuid[coreid], &cpuset);

      /* Apply this mask to the runner's pthread. */
      if (pthread_setaffinity_np(e->runners[k].thread, sizeof(cpu_set_t),
                                 &cpuset) != 0)
        error("Failed to set thread affinity.");

#else
      error("SWIFT was not compiled with affinity enabled.");
#endif
    } else {
      e->runners[k].cpuid = k;
      e->runners[k].qid = k * nr_queues / e->nr_threads;
    }

    /* Allocate particle caches. */
    e->runners[k].ci_gravity_cache.count = 0;
    e->runners[k].cj_gravity_cache.count = 0;
    gravity_cache_init(&e->runners[k].ci_gravity_cache, space_splitsize);
    gravity_cache_init(&e->runners[k].cj_gravity_cache, space_splitsize);
#ifdef WITH_VECTORIZATION
    e->runners[k].ci_cache.count = 0;
    e->runners[k].cj_cache.count = 0;
    cache_init(&e->runners[k].ci_cache, CACHE_SIZE);
    cache_init(&e->runners[k].cj_cache, CACHE_SIZE);
#endif

    if (verbose) {
      if (with_aff)
        message("runner %i on cpuid=%i with qid=%i.", e->runners[k].id,
                e->runners[k].cpuid, e->runners[k].qid);
      else
        message("runner %i using qid=%i no cpuid.", e->runners[k].id,
                e->runners[k].qid);
    }
  }

#ifdef WITH_CSDS
  if ((e->policy & engine_policy_csds) && !restart) {
    /* Write the particle csds header */
    csds_write_file_header(e->csds);
  }
#endif

  /* Initialise the structure finder */
#ifdef HAVE_VELOCIRAPTOR
  if (e->policy & engine_policy_structure_finding) velociraptor_init(e);
#endif

    /* Free the affinity stuff */
#if defined(HAVE_SETAFFINITY)
  if (with_aff) {
    free(cpuid);
  }
#endif

#ifdef SWIFT_DUMPER_THREAD

  /* Start the dumper thread.*/
  engine_dumper_init(e);
#endif

  /* Wait for the runner threads to be in place. */
  swift_barrier_wait(&e->wait_barrier);

  if (e->nodeID == 0) {
    clocks_gettime(&toc);
    message("took %.3f %s.", clocks_diff(&tic, &toc), clocks_getunit());
    fflush(stdout);
  }
}
