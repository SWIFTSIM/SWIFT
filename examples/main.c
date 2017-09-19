/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <errno.h>
#include <fenv.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Local headers. */
#include "swift.h"

/* Engine policy flags. */
#ifndef ENGINE_POLICY
#define ENGINE_POLICY engine_policy_none
#endif

/* Global profiler. */
struct profiler prof;

/**
 * @brief Help messages for the command line parameters.
 */
void print_help_message() {

  printf("\nUsage: swift [OPTION]... PARAMFILE\n");
  printf("       swift_mpi [OPTION]... PARAMFILE\n\n");

  printf("Valid options are:\n");
  printf("  %2s %14s %s\n", "-a", "", "Pin runners using processor affinity.");
  printf("  %2s %14s %s\n", "-c", "",
         "Run with cosmological time integration.");
  printf("  %2s %14s %s\n", "-C", "", "Run with cooling.");
  printf(
      "  %2s %14s %s\n", "-d", "",
      "Dry run. Read the parameter file, allocate memory but does not read ");
  printf(
      "  %2s %14s %s\n", "", "",
      "the particles from ICs and exit before the start of time integration.");
  printf("  %2s %14s %s\n", "", "",
         "Allows user to check validy of parameter and IC files as well as "
         "memory limits.");
  printf("  %2s %14s %s\n", "-D", "",
         "Always drift all particles even the ones far from active particles. "
         "This emulates");
  printf("  %2s %14s %s\n", "", "",
         "Gadget-[23] and GIZMO's default behaviours.");
  printf("  %2s %14s %s\n", "-e", "",
         "Enable floating-point exceptions (debugging mode).");
  printf("  %2s %14s %s\n", "-f", "{int}",
         "Overwrite the CPU frequency (Hz) to be used for time measurements.");
  printf("  %2s %14s %s\n", "-g", "",
         "Run with an external gravitational potential.");
  printf("  %2s %14s %s\n", "-G", "", "Run with self-gravity.");
  printf("  %2s %14s %s\n", "-M", "",
         "Reconstruct the multipoles every time-step.");
  printf("  %2s %14s %s\n", "-n", "{int}",
         "Execute a fixed number of time steps. When unset use the time_end "
         "parameter to stop.");
  printf("  %2s %14s %s\n", "-P", "{sec:par:val}",
         "Set parameter value and overwrites values read from the parameters "
         "file. Can be used more than once.");
  printf("  %2s %14s %s\n", "-s", "", "Run with hydrodynamics.");
  printf("  %2s %14s %s\n", "-S", "", "Run with stars.");
  printf("  %2s %14s %s\n", "-t", "{int}",
         "The number of threads to use on each MPI rank. Defaults to 1 if not "
         "specified.");
  printf("  %2s %14s %s\n", "-T", "", "Print timers every time-step.");
  printf("  %2s %14s %s\n", "-v", "[12]", "Increase the level of verbosity:");
  printf("  %2s %14s %s\n", "", "", "1: MPI-rank 0 writes,");
  printf("  %2s %14s %s\n", "", "", "2: All MPI-ranks write.");
  printf("  %2s %14s %s\n", "-y", "{int}",
         "Time-step frequency at which task graphs are dumped.");
  printf("  %2s %14s %s\n", "-Y", "{int}",
         "Time-step frequency at which threadpool tasks are dumped.");
  printf("  %2s %14s %s\n", "-h", "", "Print this help message and exit.");
  printf(
      "\nSee the file parameter_example.yml for an example of "
      "parameter file.\n");
}

/**
 * @brief Main routine that loads a few particles and generates some output.
 *
 */
int main(int argc, char *argv[]) {

  struct clocks_time tic, toc;

  int nr_nodes = 1, myrank = 0;

#ifdef WITH_MPI
  /* Start by initializing MPI. */
  int res = 0, prov = 0;
  if ((res = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &prov)) !=
      MPI_SUCCESS)
    error("Call to MPI_Init failed with error %i.", res);
  if (prov != MPI_THREAD_MULTIPLE)
    error(
        "MPI does not provide the level of threading"
        " required (MPI_THREAD_MULTIPLE).");
  if ((res = MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes)) != MPI_SUCCESS)
    error("MPI_Comm_size failed with error %i.", res);
  if ((res = MPI_Comm_rank(MPI_COMM_WORLD, &myrank)) != MPI_SUCCESS)
    error("Call to MPI_Comm_rank failed with error %i.", res);

  /* Make sure messages are stamped with the correct rank. */
  engine_rank = myrank;

  if ((res = MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN)) !=
      MPI_SUCCESS)
    error("Call to MPI_Comm_set_errhandler failed with error %i.", res);
  if (myrank == 0)
    printf("[0000] [00000.0] main: MPI is up and running with %i node(s).\n\n",
           nr_nodes);
  if (nr_nodes == 1) {
    message("WARNING: you are running with one MPI rank.");
    message("WARNING: you should use the non-MPI version of this program.");
  }
  fflush(stdout);

#endif

  /* Welcome to SWIFT, you made the right choice */
  if (myrank == 0) greetings();

  int with_aff = 0;
  int dry_run = 0;
  int dump_tasks = 0;
  int dump_threadpool = 0;
  int nsteps = -2;
  int with_cosmology = 0;
  int with_external_gravity = 0;
  int with_sourceterms = 0;
  int with_cooling = 0;
  int with_self_gravity = 0;
  int with_hydro = 0;
  int with_stars = 0;
  int with_fp_exceptions = 0;
  int with_drift_all = 0;
  int with_mpole_reconstruction = 0;
  int verbose = 0;
  int nr_threads = 1;
  int with_verbose_timers = 0;
  int nparams = 0;
  char *cmdparams[PARSER_MAX_NO_OF_PARAMS];
  char paramFileName[200] = "";
  unsigned long long cpufreq = 0;

  /* Parse the parameters */
  int c;
  while ((c = getopt(argc, argv, "acCdDef:FgGhMn:P:sSt:Tv:y:Y:")) != -1)
    switch (c) {
      case 'a':
#if defined(HAVE_SETAFFINITY) && defined(HAVE_LIBNUMA)
        with_aff = 1;
#else
        error("Need NUMA support for thread affinity");
#endif
        break;
      case 'c':
        with_cosmology = 1;
        break;
      case 'C':
        with_cooling = 1;
        break;
      case 'd':
        dry_run = 1;
        break;
      case 'D':
        with_drift_all = 1;
        break;
      case 'e':
        with_fp_exceptions = 1;
        break;
      case 'f':
        if (sscanf(optarg, "%llu", &cpufreq) != 1) {
          if (myrank == 0) printf("Error parsing CPU frequency (-f).\n");
          if (myrank == 0) print_help_message();
          return 1;
        }
        break;
      case 'F':
        with_sourceterms = 1;
        break;
      case 'g':
        with_external_gravity = 1;
        break;
      case 'G':
        with_self_gravity = 1;
        break;
      case 'h':
        if (myrank == 0) print_help_message();
        return 0;
      case 'M':
        with_mpole_reconstruction = 1;
        break;
      case 'n':
        if (sscanf(optarg, "%d", &nsteps) != 1) {
          if (myrank == 0) printf("Error parsing fixed number of steps.\n");
          if (myrank == 0) print_help_message();
          return 1;
        }
        break;
      case 'P':
        cmdparams[nparams] = optarg;
        nparams++;
        break;
      case 's':
        with_hydro = 1;
        break;
      case 'S':
        with_stars = 1;
        break;
      case 't':
        if (sscanf(optarg, "%d", &nr_threads) != 1) {
          if (myrank == 0)
            printf("Error parsing the number of threads (-t).\n");
          if (myrank == 0) print_help_message();
          return 1;
        }
        break;
      case 'T':
        with_verbose_timers = 1;
        break;
      case 'v':
        if (sscanf(optarg, "%d", &verbose) != 1) {
          if (myrank == 0) printf("Error parsing verbosity level (-v).\n");
          if (myrank == 0) print_help_message();
          return 1;
        }
        break;
      case 'y':
        if (sscanf(optarg, "%d", &dump_tasks) != 1) {
          if (myrank == 0) printf("Error parsing dump_tasks (-y). \n");
          if (myrank == 0) print_help_message();
          return 1;
        }
#ifndef SWIFT_DEBUG_TASKS
        if (dump_tasks) {
          error(
              "Task dumping is only possible if SWIFT was configured with the "
              "--enable-task-debugging option.");
        }
#endif
        break;
      case 'Y':
        if (sscanf(optarg, "%d", &dump_threadpool) != 1) {
          if (myrank == 0) printf("Error parsing dump_threadpool (-Y). \n");
          if (myrank == 0) print_help_message();
          return 1;
        }
#ifndef SWIFT_DEBUG_THREADPOOL
        if (dump_threadpool) {
          error(
              "Threadpool dumping is only possible if SWIFT was configured "
              "with the "
              "--enable-threadpool-debugging option.");
        }
#endif
        break;
      case '?':
        if (myrank == 0) print_help_message();
        return 1;
        break;
    }
  if (optind == argc - 1) {
    if (!strcpy(paramFileName, argv[optind++]))
      error("Error reading parameter file name.");
  } else if (optind > argc - 1) {
    if (myrank == 0) printf("Error: A parameter file name must be provided\n");
    if (myrank == 0) print_help_message();
    return 1;
  } else {
    if (myrank == 0) printf("Error: Too many parameters given\n");
    if (myrank == 0) print_help_message();
    return 1;
  }
  if (!with_self_gravity && !with_hydro && !with_external_gravity) {
    if (myrank == 0)
      printf("Error: At least one of -s, -g or -G must be chosen.\n");
    if (myrank == 0) print_help_message();
    return 1;
  }
  if (with_stars && !with_external_gravity && !with_self_gravity) {
    if (myrank == 0)
      printf(
          "Error: Cannot process stars without gravity, -g or -G must be "
          "chosen.\n");
    if (myrank == 0) print_help_message();
    return 1;
  }

/* Let's pin the main thread, now we know if affinity will be used. */
#if defined(HAVE_SETAFFINITY) && defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
  if (with_aff && ((ENGINE_POLICY)&engine_policy_setaffinity) == engine_policy_setaffinity)
    engine_pin();
#endif

  /* Genesis 1.1: And then, there was time ! */
  clocks_set_cpufreq(cpufreq);

  /* How vocal are we ? */
  const int talking = (verbose == 1 && myrank == 0) || (verbose == 2);

  if (myrank == 0 && dry_run)
    message(
        "Executing a dry run. No i/o or time integration will be performed.");

  /* Report CPU frequency.*/
  cpufreq = clocks_get_cpufreq();
  if (myrank == 0) {
    message("CPU frequency used for tick conversion: %llu Hz", cpufreq);
  }

/* Report host name(s). */
#ifdef WITH_MPI
  if (talking) {
    message("Rank %d running on: %s", myrank, hostname());
  }
#else
  message("Running on: %s", hostname());
#endif

/* Do we have debugging checks ? */
#ifdef SWIFT_DEBUG_CHECKS
  if (myrank == 0)
    message("WARNING: Debugging checks activated. Code will be slower !");
#endif

/* Do we have gravity accuracy checks ? */
#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  if (myrank == 0)
    message(
        "WARNING: Checking 1/%d of all gpart for gravity accuracy. Code will "
        "be slower !",
        SWIFT_GRAVITY_FORCE_CHECKS);
#endif

  /* Do we choke on FP-exceptions ? */
  if (with_fp_exceptions) {
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    if (myrank == 0)
      message("WARNING: Floating point exceptions will be reported.");
  }

  /* How large are the parts? */
  if (myrank == 0) {
    message("sizeof(part)        is %4zi bytes.", sizeof(struct part));
    message("sizeof(xpart)       is %4zi bytes.", sizeof(struct xpart));
    message("sizeof(spart)       is %4zi bytes.", sizeof(struct spart));
    message("sizeof(gpart)       is %4zi bytes.", sizeof(struct gpart));
    message("sizeof(multipole)   is %4zi bytes.", sizeof(struct multipole));
    message("sizeof(grav_tensor) is %4zi bytes.", sizeof(struct grav_tensor));
    message("sizeof(task)        is %4zi bytes.", sizeof(struct task));
    message("sizeof(cell)        is %4zi bytes.", sizeof(struct cell));
  }

  /* Read the parameter file */
  struct swift_params *params = malloc(sizeof(struct swift_params));
  if (params == NULL) error("Error allocating memory for the parameter file.");
  if (myrank == 0) {
    message("Reading runtime parameters from file '%s'", paramFileName);
    parser_read_file(paramFileName, params);

    /* Handle any command-line overrides. */
    if (nparams > 0) {
      message(
          "Overwriting values read from the YAML file with command-line "
          "values.");
      for (int k = 0; k < nparams; k++) parser_set_param(params, cmdparams[k]);
    }

    /* And dump the parameters as used. */
    // parser_print_params(&params);
    parser_write_params_to_file(params, "used_parameters.yml");
  }
#ifdef WITH_MPI
  /* Broadcast the parameter file */
  MPI_Bcast(params, sizeof(struct swift_params), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  /* Check that we can write the snapshots by testing if the output
   * directory exists and is searchable and writable. */
  char basename[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "Snapshots:basename", basename);
  const char *dirp = dirname(basename);
  if (access(dirp, W_OK | X_OK) != 0) {
    error("Cannot write snapshots in directory %s (%s)", dirp, strerror(errno));
  }

  /* Prepare the domain decomposition scheme */
  struct repartition reparttype;
#ifdef WITH_MPI
  struct partition initial_partition;
  partition_init(&initial_partition, &reparttype, params, nr_nodes);

  /* Let's report what we did */
  if (myrank == 0) {
    message("Using initial partition %s",
            initial_partition_name[initial_partition.type]);
    if (initial_partition.type == INITPART_GRID)
      message("grid set to [ %i %i %i ].", initial_partition.grid[0],
              initial_partition.grid[1], initial_partition.grid[2]);
    message("Using %s repartitioning", repartition_name[reparttype.type]);
  }
#endif

  /* Initialize unit system and constants */
  struct unit_system us;
  struct phys_const prog_const;
  units_init(&us, params, "InternalUnitSystem");
  phys_const_init(&us, &prog_const);
  if (myrank == 0 && verbose > 0) {
    message("Internal unit system: U_M = %e g.", us.UnitMass_in_cgs);
    message("Internal unit system: U_L = %e cm.", us.UnitLength_in_cgs);
    message("Internal unit system: U_t = %e s.", us.UnitTime_in_cgs);
    message("Internal unit system: U_I = %e A.", us.UnitCurrent_in_cgs);
    message("Internal unit system: U_T = %e K.", us.UnitTemperature_in_cgs);
    phys_const_print(&prog_const);
  }

  /* Initialise the hydro properties */
  struct hydro_props hydro_properties;
  if (with_hydro) hydro_props_init(&hydro_properties, params);

  /* Initialise the gravity properties */
  struct gravity_props gravity_properties;
  if (with_self_gravity) gravity_props_init(&gravity_properties, params);

  /* Read particles and space information from (GADGET) ICs */
  char ICfileName[200] = "";
  parser_get_param_string(params, "InitialConditions:file_name", ICfileName);
  const int replicate =
      parser_get_opt_param_int(params, "InitialConditions:replicate", 1);
  const int clean_h_values =
      parser_get_opt_param_int(params, "InitialConditions:cleanup_h", 0);
  if (myrank == 0) message("Reading ICs from file '%s'", ICfileName);
  fflush(stdout);

  /* Get ready to read particles of all kinds */
  struct part *parts = NULL;
  struct gpart *gparts = NULL;
  struct spart *sparts = NULL;
  size_t Ngas = 0, Ngpart = 0, Nspart = 0;
  double dim[3] = {0., 0., 0.};
  int periodic = 0;
  int flag_entropy_ICs = 0;
  if (myrank == 0) clocks_gettime(&tic);
#if defined(WITH_MPI)
#if defined(HAVE_PARALLEL_HDF5)
  read_ic_parallel(ICfileName, &us, dim, &parts, &gparts, &sparts, &Ngas,
                   &Ngpart, &Nspart, &periodic, &flag_entropy_ICs, with_hydro,
                   (with_external_gravity || with_self_gravity), with_stars,
                   myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL, dry_run);
#else
  read_ic_serial(ICfileName, &us, dim, &parts, &gparts, &sparts, &Ngas, &Ngpart,
                 &Nspart, &periodic, &flag_entropy_ICs, with_hydro,
                 (with_external_gravity || with_self_gravity), with_stars,
                 myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL, dry_run);
#endif
#else
  read_ic_single(ICfileName, &us, dim, &parts, &gparts, &sparts, &Ngas, &Ngpart,
                 &Nspart, &periodic, &flag_entropy_ICs, with_hydro,
                 (with_external_gravity || with_self_gravity), with_stars,
                 dry_run);
#endif
  if (myrank == 0) {
    clocks_gettime(&toc);
    message("Reading initial conditions took %.3f %s.", clocks_diff(&tic, &toc),
            clocks_getunit());
    fflush(stdout);
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check once and for all that we don't have unwanted links */
  if (!with_stars) {
    for (size_t k = 0; k < Ngpart; ++k)
      if (gparts[k].type == swift_type_star) error("Linking problem");
  }
  if (!with_hydro) {
    for (size_t k = 0; k < Ngpart; ++k)
      if (gparts[k].type == swift_type_gas) error("Linking problem");
  }
#endif

  /* Get the total number of particles across all nodes. */
  long long N_total[3] = {0, 0, 0};
#if defined(WITH_MPI)
  long long N_long[3] = {Ngas, Ngpart, Nspart};
  MPI_Allreduce(&N_long, &N_total, 3, MPI_LONG_LONG_INT, MPI_SUM,
                MPI_COMM_WORLD);
#else
  N_total[0] = Ngas;
  N_total[1] = Ngpart;
  N_total[2] = Nspart;
#endif
  if (myrank == 0)
    message(
        "Read %lld gas particles, %lld star particles and %lld gparts from the "
        "ICs.",
        N_total[0], N_total[2], N_total[1]);

  /* Initialize the space with these data. */
  if (myrank == 0) clocks_gettime(&tic);
  struct space s;
  space_init(&s, params, dim, parts, gparts, sparts, Ngas, Ngpart, Nspart,
             periodic, replicate, with_self_gravity, talking, dry_run);
  if (myrank == 0) {
    clocks_gettime(&toc);
    message("space_init took %.3f %s.", clocks_diff(&tic, &toc),
            clocks_getunit());
    fflush(stdout);
  }

  /* Say a few nice things about the space we just created. */
  if (myrank == 0) {
    message("space dimensions are [ %.3f %.3f %.3f ].", s.dim[0], s.dim[1],
            s.dim[2]);
    message("space %s periodic.", s.periodic ? "is" : "isn't");
    message("highest-level cell dimensions are [ %i %i %i ].", s.cdim[0],
            s.cdim[1], s.cdim[2]);
    message("%zi parts in %i cells.", s.nr_parts, s.tot_cells);
    message("%zi gparts in %i cells.", s.nr_gparts, s.tot_cells);
    message("%zi sparts in %i cells.", s.nr_sparts, s.tot_cells);
    message("maximum depth is %d.", s.maxdepth);
    fflush(stdout);
  }

  /* Verify that each particle is in it's proper cell. */
  if (talking && !dry_run) {
    int icount = 0;
    space_map_cells_pre(&s, 0, &map_cellcheck, &icount);
    message("map_cellcheck picked up %i parts.", icount);
  }

  /* Verify the maximal depth of cells. */
  if (talking && !dry_run) {
    int data[2] = {s.maxdepth, 0};
    space_map_cells_pre(&s, 0, &map_maxdepth, data);
    message("nr of cells at depth %i is %i.", data[0], data[1]);
  }

/* Initialise the table of Ewald corrections for the gravity checks */
#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  if (periodic) gravity_exact_force_ewald_init(dim[0]);
#endif

  /* Initialise the external potential properties */
  struct external_potential potential;
  if (with_external_gravity)
    potential_init(params, &prog_const, &us, &s, &potential);
  if (with_external_gravity && myrank == 0) potential_print(&potential);

  /* Initialise the cooling function properties */
  struct cooling_function_data cooling_func;
  if (with_cooling) cooling_init(params, &us, &prog_const, &cooling_func);
  if (with_cooling && myrank == 0) cooling_print(&cooling_func);

  /* Initialise the feedback properties */
  struct sourceterms sourceterms;
  if (with_sourceterms) sourceterms_init(params, &us, &sourceterms);
  if (with_sourceterms && myrank == 0) sourceterms_print(&sourceterms);

  /* Construct the engine policy */
  int engine_policies = ENGINE_POLICY | engine_policy_steal;
  if (with_drift_all) engine_policies |= engine_policy_drift_all;
  if (with_mpole_reconstruction)
    engine_policies |= engine_policy_reconstruct_mpoles;
  if (with_hydro) engine_policies |= engine_policy_hydro;
  if (with_self_gravity) engine_policies |= engine_policy_self_gravity;
  if (with_external_gravity) engine_policies |= engine_policy_external_gravity;
  if (with_cosmology) engine_policies |= engine_policy_cosmology;
  if (with_cooling) engine_policies |= engine_policy_cooling;
  if (with_sourceterms) engine_policies |= engine_policy_sourceterms;
  if (with_stars) engine_policies |= engine_policy_stars;

  /* Initialize the engine with the space and policies. */
  if (myrank == 0) clocks_gettime(&tic);
  struct engine e;
  engine_init(&e, &s, params, nr_nodes, myrank, nr_threads, N_total[0],
              N_total[1], with_aff, engine_policies, talking, &reparttype, &us,
              &prog_const, &hydro_properties, &gravity_properties, &potential,
              &cooling_func, &sourceterms);
  if (myrank == 0) {
    clocks_gettime(&toc);
    message("engine_init took %.3f %s.", clocks_diff(&tic, &toc),
            clocks_getunit());
    fflush(stdout);
  }

/* Init the runner history. */
#ifdef HIST
  for (k = 0; k < runner_hist_N; k++) runner_hist_bins[k] = 0;
#endif

#if defined(WITH_MPI)
  N_long[0] = s.nr_parts;
  N_long[1] = s.nr_gparts;
  N_long[2] = s.nr_sparts;
  MPI_Reduce(&N_long, &N_total, 3, MPI_LONG_LONG_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);
#else
  N_total[0] = s.nr_parts;
  N_total[1] = s.nr_gparts;
  N_total[2] = s.nr_sparts;
#endif

  /* Get some info to the user. */
  if (myrank == 0) {
    long long N_DM = N_total[1] - N_total[2] - N_total[0];
    message(
        "Running on %lld gas particles, %lld star particles and %lld DM "
        "particles (%lld gravity particles)",
        N_total[0], N_total[2], N_total[1] > 0 ? N_DM : 0, N_total[1]);
    message(
        "from t=%.3e until t=%.3e with %d threads and %d queues (dt_min=%.3e, "
        "dt_max=%.3e)...",
        e.timeBegin, e.timeEnd, e.nr_threads, e.sched.nr_queues, e.dt_min,
        e.dt_max);
    fflush(stdout);
  }

  /* Time to say good-bye if this was not a serious run. */
  if (dry_run) {
#ifdef WITH_MPI
    if ((res = MPI_Finalize()) != MPI_SUCCESS)
      error("call to MPI_Finalize failed with error %i.", res);
#endif
    if (myrank == 0)
      message("Time integration ready to start. End of dry-run.");
    engine_clean(&e);
    free(params);
    return 0;
  }

#ifdef WITH_MPI
  /* Split the space. */
  engine_split(&e, &initial_partition);
  engine_redistribute(&e);
#endif

  /* Initialise the particles */
  engine_init_particles(&e, flag_entropy_ICs, clean_h_values);

  /* Write the state of the system before starting time integration. */
  engine_dump_snapshot(&e);
  engine_print_stats(&e);

  /* Legend */
  if (myrank == 0)
    printf("# %6s %14s %14s %10s %10s %10s %16s [%s]\n", "Step", "Time",
           "Time-step", "Updates", "g-Updates", "s-Updates", "Wall-clock time",
           clocks_getunit());

  /* File for the timers */
  if (with_verbose_timers) timers_open_file(myrank);

  /* Main simulation loop */
  for (int j = 0; !engine_is_done(&e) && e.step - 1 != nsteps; j++) {

    /* Reset timers */
    timers_reset_all();

    /* Take a step. */
    engine_step(&e);

    /* Print the timers. */
    if (with_verbose_timers) timers_print(e.step);

#ifdef SWIFT_DEBUG_TASKS
    /* Dump the task data using the given frequency. */
    if (dump_tasks && (dump_tasks == 1 || j % dump_tasks == 1)) {
#ifdef WITH_MPI

      /* Make sure output file is empty, only on one rank. */
      char dumpfile[30];
      snprintf(dumpfile, 30, "thread_info_MPI-step%d.dat", j + 1);
      FILE *file_thread;
      if (myrank == 0) {
        file_thread = fopen(dumpfile, "w");
        fclose(file_thread);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      for (int i = 0; i < nr_nodes; i++) {

        /* Rank 0 decides the index of writing node, this happens one-by-one. */
        int kk = i;
        MPI_Bcast(&kk, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (i == myrank) {

          /* Open file and position at end. */
          file_thread = fopen(dumpfile, "a");

          fprintf(file_thread, " %03i 0 0 0 0 %lli %lli %zi %zi %zi 0 0 %lli\n",
                  myrank, e.tic_step, e.toc_step, e.updates, e.g_updates,
                  e.s_updates, cpufreq);
          int count = 0;
          for (int l = 0; l < e.sched.nr_tasks; l++) {
            if (!e.sched.tasks[l].implicit && e.sched.tasks[l].toc != 0) {
              fprintf(
                  file_thread,
                  " %03i %i %i %i %i %lli %lli %i %i %i %i %i %i\n", myrank,
                  e.sched.tasks[l].rid, e.sched.tasks[l].type,
                  e.sched.tasks[l].subtype, (e.sched.tasks[l].cj == NULL),
                  e.sched.tasks[l].tic, e.sched.tasks[l].toc,
                  (e.sched.tasks[l].ci != NULL) ? e.sched.tasks[l].ci->count
                                                : 0,
                  (e.sched.tasks[l].cj != NULL) ? e.sched.tasks[l].cj->count
                                                : 0,
                  (e.sched.tasks[l].ci != NULL) ? e.sched.tasks[l].ci->gcount
                                                : 0,
                  (e.sched.tasks[l].cj != NULL) ? e.sched.tasks[l].cj->gcount
                                                : 0,
                  e.sched.tasks[l].flags, e.sched.tasks[l].sid);
            }
            fflush(stdout);
            count++;
          }
          fclose(file_thread);
        }

        /* And we wait for all to synchronize. */
        MPI_Barrier(MPI_COMM_WORLD);
      }

#else
      char dumpfile[30];
      snprintf(dumpfile, 30, "thread_info-step%d.dat", j + 1);
      FILE *file_thread;
      file_thread = fopen(dumpfile, "w");
      /* Add some information to help with the plots */
      fprintf(file_thread, " %i %i %i %i %lli %lli %zi %zi %zi %i %lli\n", -2,
              -1, -1, 1, e.tic_step, e.toc_step, e.updates, e.g_updates,
              e.s_updates, 0, cpufreq);
      for (int l = 0; l < e.sched.nr_tasks; l++) {
        if (!e.sched.tasks[l].implicit && e.sched.tasks[l].toc != 0) {
          fprintf(
              file_thread, " %i %i %i %i %lli %lli %i %i %i %i %i\n",
              e.sched.tasks[l].rid, e.sched.tasks[l].type,
              e.sched.tasks[l].subtype, (e.sched.tasks[l].cj == NULL),
              e.sched.tasks[l].tic, e.sched.tasks[l].toc,
              (e.sched.tasks[l].ci == NULL) ? 0 : e.sched.tasks[l].ci->count,
              (e.sched.tasks[l].cj == NULL) ? 0 : e.sched.tasks[l].cj->count,
              (e.sched.tasks[l].ci == NULL) ? 0 : e.sched.tasks[l].ci->gcount,
              (e.sched.tasks[l].cj == NULL) ? 0 : e.sched.tasks[l].cj->gcount,
              e.sched.tasks[l].sid);
        }
      }
      fclose(file_thread);
#endif  // WITH_MPI
    }
#endif  // SWIFT_DEBUG_TASKS

#ifdef SWIFT_DEBUG_THREADPOOL
    /* Dump the task data using the given frequency. */
    if (dump_threadpool && (dump_threadpool == 1 || j % dump_threadpool == 1)) {
      char dumpfile[40];
#ifdef WITH_MPI
      snprintf(dumpfile, 30, "threadpool_info-rank%d-step%d.dat", engine_rank,
               j + 1);
#else
      snprintf(dumpfile, 30, "threadpool_info-step%d.dat", j + 1);
#endif  // WITH_MPI
      threadpool_dump_log(&e.threadpool, dumpfile, 1);
    } else {
      threadpool_reset_log(&e.threadpool);
    }
#endif  // SWIFT_DEBUG_THREADPOOL
  }

/* Print the values of the runner histogram. */
#ifdef HIST
  printf("main: runner histogram data:\n");
  for (k = 0; k < runner_hist_N; k++)
    printf(" %e %e %e\n",
           runner_hist_a + k * (runner_hist_b - runner_hist_a) / runner_hist_N,
           runner_hist_a +
               (k + 1) * (runner_hist_b - runner_hist_a) / runner_hist_N,
           (double)runner_hist_bins[k]);
#endif

  /* Write final output. */
  engine_drift_all(&e);
  engine_print_stats(&e);
  engine_dump_snapshot(&e);

#ifdef WITH_MPI
  if ((res = MPI_Finalize()) != MPI_SUCCESS)
    error("call to MPI_Finalize failed with error %i.", res);
#endif

  /* Clean everything */
  if (with_verbose_timers) timers_close_file();
  engine_clean(&e);
  free(params);

  /* Say goodbye. */
  if (myrank == 0) message("done. Bye.");

  /* All is calm, all is bright. */
  return 0;
}
