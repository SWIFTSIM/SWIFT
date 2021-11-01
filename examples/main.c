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
#include <sys/stat.h>
#include <unistd.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* Local headers. */
#include "argparse.h"
#include "swift.h"

/* Engine policy flags. */
#ifndef ENGINE_POLICY
#define ENGINE_POLICY engine_policy_none
#endif

/* Global profiler. */
struct profiler prof;

/*  Usage string. */
static const char *const swift_usage[] = {
    "swift [options] [[--] param-file]",
    "swift [options] param-file",
    "swift_mpi [options] [[--] param-file]",
    "swift_mpi [options] param-file",
    NULL,
};

/* Function to handle multiple -P arguments. */
struct cmdparams {
  const char *param[PARSER_MAX_NO_OF_PARAMS];
  int nparam;
};

static int handle_cmdparam(struct argparse *self,
                           const struct argparse_option *opt) {
  struct cmdparams *cmdps = (struct cmdparams *)opt->data;
  cmdps->param[cmdps->nparam] = *(char **)opt->value;
  cmdps->nparam++;
  return 1;
}

/**
 * @brief Main routine that loads a few particles and generates some output.
 *
 */
int main(int argc, char *argv[]) {

  struct clocks_time tic, toc;
  struct engine e;

  /* Structs used by the engine. Declare now to make sure these are always in
   * scope.  */
  struct chemistry_global_data chemistry;
  struct cooling_function_data cooling_func;
  struct cosmology cosmo;
  struct external_potential potential;
  struct star_formation starform;
  struct pm_mesh mesh;
  struct gpart *gparts = NULL;
  struct gravity_props gravity_properties;
  struct hydro_props hydro_properties;
  struct stars_props stars_properties;
  struct sink_props sink_properties;
  struct neutrino_props neutrino_properties;
  struct feedback_props feedback_properties;
  struct rt_props rt_properties;
  struct entropy_floor_properties entropy_floor;
  struct black_holes_props black_holes_properties;
  struct fof_props fof_properties;
  struct part *parts = NULL;
  struct phys_const prog_const;
  struct space s;
  struct spart *sparts = NULL;
  struct bpart *bparts = NULL;
  struct sink *sinks = NULL;
  struct unit_system us;
  struct los_props los_properties;
  struct ic_info ics_metadata;

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

  /* Make sure messages are stamped with the correct rank and step. */
  engine_rank = myrank;
  engine_current_step = 0;

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
  if (myrank == 0) greetings(/*fof=*/0);

  int with_aff = 0;
  int dry_run = 0;
  int dump_tasks = 0;
  int dump_cells = 0;
  int dump_threadpool = 0;
  int nsteps = -2;
  int restart = 0;
  int with_cosmology = 0;
  int with_external_gravity = 0;
  int with_temperature = 0;
  int with_cooling = 0;
  int with_self_gravity = 0;
  int with_hydro = 0;
  int with_stars = 0;
  int with_fof = 0;
  int with_star_formation = 0;
  int with_feedback = 0;
  int with_black_holes = 0;
  int with_timestep_limiter = 0;
  int with_timestep_sync = 0;
  int with_fp_exceptions = 0;
  int with_drift_all = 0;
  int with_mpole_reconstruction = 0;
  int with_structure_finding = 0;
  int with_csds = 0;
  int with_sink = 0;
  int with_qla = 0;
  int with_eagle = 0;
  int with_gear = 0;
  int with_line_of_sight = 0;
  int with_rt = 0;
  int verbose = 0;
  int nr_threads = 1;
  int nr_pool_threads = -1;
  int with_verbose_timers = 0;
  char *output_parameters_filename = NULL;
  char *cpufreqarg = NULL;
  char *param_filename = NULL;
  char restart_file[200] = "";
  unsigned long long cpufreq = 0;
  float dump_tasks_threshold = 0.f;
  struct cmdparams cmdps;
  cmdps.nparam = 0;
  cmdps.param[0] = NULL;
  char *buffer = NULL;

  /* Parse the command-line parameters. */
  struct argparse_option options[] = {
      OPT_HELP(),

      OPT_GROUP("  Simulation options:\n"),
      OPT_BOOLEAN('b', "feedback", &with_feedback, "Run with stars feedback.",
                  NULL, 0, 0),
      OPT_BOOLEAN('c', "cosmology", &with_cosmology,
                  "Run with cosmological time integration.", NULL, 0, 0),
      OPT_BOOLEAN(0, "temperature", &with_temperature,
                  "Run with temperature calculation.", NULL, 0, 0),
      OPT_BOOLEAN('C', "cooling", &with_cooling,
                  "Run with cooling (also switches on --temperature).", NULL, 0,
                  0),
      OPT_BOOLEAN('D', "drift-all", &with_drift_all,
                  "Always drift all particles even the ones far from active "
                  "particles. This emulates Gadget-[23] and GIZMO's default "
                  "behaviours.",
                  NULL, 0, 0),
      OPT_BOOLEAN('F', "star-formation", &with_star_formation,
                  "Run with star formation.", NULL, 0, 0),
      OPT_BOOLEAN('g', "external-gravity", &with_external_gravity,
                  "Run with an external gravitational potential.", NULL, 0, 0),
      OPT_BOOLEAN('G', "self-gravity", &with_self_gravity,
                  "Run with self-gravity.", NULL, 0, 0),
      OPT_BOOLEAN('M', "multipole-reconstruction", &with_mpole_reconstruction,
                  "Reconstruct the multipoles every time-step.", NULL, 0, 0),
      OPT_BOOLEAN('s', "hydro", &with_hydro, "Run with hydrodynamics.", NULL, 0,
                  0),
      OPT_BOOLEAN('S', "stars", &with_stars, "Run with stars.", NULL, 0, 0),
      OPT_BOOLEAN('B', "black-holes", &with_black_holes,
                  "Run with black holes.", NULL, 0, 0),
      OPT_BOOLEAN('k', "sinks", &with_sink, "Run with sink particles.", NULL, 0,
                  0),
      OPT_BOOLEAN(
          'u', "fof", &with_fof,
          "Run Friends-of-Friends algorithm to perform black hole seeding.",
          NULL, 0, 0),
      OPT_BOOLEAN('x', "velociraptor", &with_structure_finding,
                  "Run with structure finding.", NULL, 0, 0),
      OPT_BOOLEAN(0, "line-of-sight", &with_line_of_sight,
                  "Run with line-of-sight outputs.", NULL, 0, 0),
      OPT_BOOLEAN(0, "limiter", &with_timestep_limiter,
                  "Run with time-step limiter.", NULL, 0, 0),
      OPT_BOOLEAN(0, "sync", &with_timestep_sync,
                  "Run with time-step synchronization of particles hit by "
                  "feedback events.",
                  NULL, 0, 0),
      OPT_BOOLEAN(0, "csds", &with_csds,
                  "Run with the Continuous Simulation Data Stream (CSDS).",
                  NULL, 0, 0),
      OPT_BOOLEAN('R', "radiation", &with_rt,
                  "Run with radiative transfer. Work in progress, currently "
                  "has no effect.",
                  NULL, 0, 0),

      OPT_GROUP("  Simulation meta-options:\n"),
      OPT_BOOLEAN(0, "quick-lyman-alpha", &with_qla,
                  "Run with all the options needed for the quick Lyman-alpha "
                  "model. This is equivalent to --hydro --self-gravity --stars "
                  "--star-formation --cooling.",
                  NULL, 0, 0),
      OPT_BOOLEAN(
          0, "eagle", &with_eagle,
          "Run with all the options needed for the EAGLE model. This is "
          "equivalent to --hydro --limiter --sync --self-gravity --stars "
          "--star-formation --cooling --feedback --black-holes --fof.",
          NULL, 0, 0),
      OPT_BOOLEAN(
          0, "gear", &with_gear,
          "Run with all the options needed for the GEAR model. This is "
          "equivalent to --hydro --limiter --sync --self-gravity --stars "
          "--star-formation --cooling --feedback.",
          NULL, 0, 0),

      OPT_GROUP("  Control options:\n"),
      OPT_BOOLEAN('a', "pin", &with_aff,
                  "Pin runners using processor affinity.", NULL, 0, 0),
      OPT_BOOLEAN('d', "dry-run", &dry_run,
                  "Dry run. Read the parameter file, allocates memory but does "
                  "not read the particles from ICs. Exits before the start of "
                  "time integration. Checks the validity of parameters and IC "
                  "files as well as memory limits.",
                  NULL, 0, 0),
      OPT_BOOLEAN('e', "fpe", &with_fp_exceptions,
                  "Enable floating-point exceptions (debugging mode).", NULL, 0,
                  0),
      OPT_STRING('f', "cpu-frequency", &cpufreqarg,
                 "Overwrite the CPU "
                 "frequency (Hz) to be used for time measurements.",
                 NULL, 0, 0),
      OPT_INTEGER('n', "steps", &nsteps,
                  "Execute a fixed number of time steps. When unset use the "
                  "time_end parameter to stop.",
                  NULL, 0, 0),
      OPT_STRING('o', "output-params", &output_parameters_filename,
                 "Generate a parameter file with the options for selecting the "
                 "output fields.",
                 NULL, 0, 0),
      OPT_STRING('P', "param", &buffer,
                 "Set parameter value, overiding the value read from the "
                 "parameter file. Can be used more than once {sec:par:value}.",
                 handle_cmdparam, (intptr_t)&cmdps, 0),
      OPT_BOOLEAN('r', "restart", &restart, "Continue using restart files.",
                  NULL, 0, 0),
      OPT_INTEGER(
          't', "threads", &nr_threads,
          "The number of task threads to use on each MPI rank. Defaults to "
          "1 if not specified.",
          NULL, 0, 0),
      OPT_INTEGER(0, "pool-threads", &nr_pool_threads,
                  "The number of threads to use on each MPI rank for the "
                  "threadpool operations. "
                  "Defaults to the numbers of task threads if not specified.",
                  NULL, 0, 0),
      OPT_INTEGER('T', "timers", &with_verbose_timers,
                  "Print timers every time-step.", NULL, 0, 0),
      OPT_INTEGER('v', "verbose", &verbose,
                  "Run in verbose mode, in MPI mode 2 outputs from all ranks.",
                  NULL, 0, 0),
      OPT_INTEGER('y', "task-dumps", &dump_tasks,
                  "Time-step frequency at which task graphs are dumped.", NULL,
                  0, 0),
      OPT_INTEGER(0, "cell-dumps", &dump_cells,
                  "Time-step frequency at which cell graphs are dumped.", NULL,
                  0, 0),
      OPT_INTEGER('Y', "threadpool-dumps", &dump_threadpool,
                  "Time-step frequency at which threadpool tasks are dumped.",
                  NULL, 0, 0),
      OPT_FLOAT(0, "dump-tasks-threshold", &dump_tasks_threshold,
                "Fraction of the total step's time spent in a task to trigger "
                "a dump of the task plot on this step",
                NULL, 0, 0),
      OPT_END(),
  };
  struct argparse argparse;
  argparse_init(&argparse, options, swift_usage, 0);
  argparse_describe(&argparse, "\nParameters:",
                    "\nSee the file examples/parameter_example.yml for an "
                    "example of parameter file.");
  int nargs = argparse_parse(&argparse, argc, (const char **)argv);

  /* Deal with meta options */
  if (with_qla) {
    with_hydro = 1;
    with_self_gravity = 1;
    with_stars = 1;
    with_star_formation = 1;
    with_cooling = 1;
  }
  if (with_eagle) {
    with_hydro = 1;
    with_timestep_limiter = 1;
    with_timestep_sync = 1;
    with_self_gravity = 1;
    with_stars = 1;
    with_star_formation = 1;
    with_cooling = 1;
    with_feedback = 1;
    with_black_holes = 1;
    with_fof = 1;
  }
  if (with_gear) {
    with_hydro = 1;
    with_timestep_limiter = 1;
    with_timestep_sync = 1;
    with_self_gravity = 1;
    with_stars = 1;
    with_star_formation = 1;
    with_cooling = 1;
    with_feedback = 1;
  }

  /* Deal with thread numbers */
  if (nr_pool_threads == -1) nr_pool_threads = nr_threads;

  /* Write output parameter file */
  if (myrank == 0 && output_parameters_filename != NULL) {
    io_write_output_field_parameter(output_parameters_filename, with_cosmology,
                                    with_fof, with_structure_finding);
    printf("End of run.\n");
    return 0;
  }

  /* Need a parameter file. */
  if (nargs != 1) {
    if (myrank == 0) argparse_usage(&argparse);
    printf("\nError: no parameter file was supplied.\n");
    return 1;
  }
  param_filename = argv[0];

  /* Checks of options. */
#if !defined(HAVE_SETAFFINITY) || !defined(HAVE_LIBNUMA)
  if (with_aff) {
    printf("Error: no NUMA support for thread affinity\n");
    return 1;
  }
#endif

#if !defined(WITH_CSDS)
  if (with_csds) {
    printf(
        "Error: the CSDS is not available, please compile with "
        "--enable-csds.");
    return 1;
  }
#endif

#ifndef HAVE_FE_ENABLE_EXCEPT
  if (with_fp_exceptions) {
    printf("Error: no support for floating point exceptions\n");
    return 1;
  }
#endif

#ifndef HAVE_VELOCIRAPTOR
  if (with_structure_finding) {
    printf("Error: VELOCIraptor is not available\n");
    return 1;
  }
#endif

#ifndef SWIFT_DEBUG_TASKS
  if (dump_tasks) {
    if (myrank == 0) {
      message(
          "WARNING: complete task dumps are only created when "
          "configured with --enable-task-debugging.");
      message("         Basic task statistics will be output.");
    }
  }
#endif

  if (dump_tasks_threshold > 0.f) {
#ifndef SWIFT_DEBUG_TASKS
    if (myrank == 0) {
      error(
          "Error: Dumping task plot data above a fixed time threshold is only "
          "valid when the code is configured with --enable-task-debugging.");
    }
#endif
#ifdef WITH_MPI
    if (nr_nodes > 1)
      error("Cannot dump tasks above a time threshold over MPI (yet).");
#endif
  }

#ifndef SWIFT_CELL_GRAPH
  if (dump_cells) {
    if (myrank == 0) {
      error(
          "complete cell dumps are only created when "
          "configured with --enable-cell-graph.");
    }
  }
#endif

#ifndef SWIFT_DEBUG_THREADPOOL
  if (dump_threadpool) {
    printf(
        "Error: threadpool dumping is only possible if SWIFT was "
        "configured with the --enable-threadpool-debugging option.\n");
    return 1;
  }
#endif

#ifdef WITH_MPI
  if (with_sink) {
    printf("Error: sink particles are not available yet with MPI.\n");
    return 1;
  }
#endif

  if (with_sink && with_star_formation) {
    printf(
        "Error: The sink particles are not supposed to be run with star "
        "formation.\n");
    return 1;
  }

  /* The CPU frequency is a long long, so we need to parse that ourselves. */
  if (cpufreqarg != NULL) {
    if (sscanf(cpufreqarg, "%llu", &cpufreq) != 1) {
      if (myrank == 0)
        printf("Error parsing CPU frequency (%s).\n", cpufreqarg);
      return 1;
    }
  }

  if (!with_self_gravity && !with_hydro && !with_external_gravity) {
    if (myrank == 0) {
      argparse_usage(&argparse);
      printf("\nError: At least one of -s, -g or -G must be chosen.\n");
    }
    return 1;
  }
  if (with_stars && !with_external_gravity && !with_self_gravity) {
    if (myrank == 0) {
      argparse_usage(&argparse);
      printf(
          "\nError: Cannot process stars without gravity, -g or -G "
          "must be chosen.\n");
    }
    return 1;
  }

  if (with_black_holes && !with_self_gravity) {
    if (myrank == 0) {
      argparse_usage(&argparse);
      printf(
          "\nError: Cannot process black holes without self-gravity, -G must "
          "be chosen.\n");
    }
    return 1;
  }

  if (with_fof) {
#ifndef WITH_FOF
    error("Running with FOF but compiled without it!");
#endif
  }

  if (with_fof && !with_self_gravity) {
    if (myrank == 0)
      printf(
          "Error: Cannot perform FOF search without gravity, -g or -G must be "
          "chosen.\n");
    return 1;
  }

  if (!with_stars && with_star_formation) {
    if (myrank == 0) {
      argparse_usage(&argparse);
      printf(
          "\nError: Cannot process star formation without stars, --stars must "
          "be chosen.\n");
    }
    return 1;
  }

  if (!with_stars && with_feedback) {
    if (myrank == 0) {
      argparse_usage(&argparse);
      printf(
          "\nError: Cannot process feedback without stars, --stars must be "
          "chosen.\n");
    }
    return 1;
  }

  if (!with_hydro && with_feedback) {
    if (myrank == 0) {
      argparse_usage(&argparse);
      printf(
          "\nError: Cannot process feedback without gas, --hydro must be "
          "chosen.\n");
    }
    return 1;
  }

  if (!with_hydro && with_black_holes) {
    if (myrank == 0) {
      argparse_usage(&argparse);
      printf(
          "\nError: Cannot process black holes without gas, --hydro must be "
          "chosen.\n");
    }
    return 1;
  }

  if (!with_hydro && with_line_of_sight) {
    if (myrank == 0) {
      argparse_usage(&argparse);
      printf(
          "\nError: Cannot use line-of-sight outputs without gas, --hydro must "
          "be chosen.\n");
    }
    return 1;
  }

#ifdef RT_NONE
  if (with_rt) {
    error("Running with radiative transfer but compiled without it!");
  }
#else
  if (with_rt && !with_hydro) {
    error(
        "Error: Cannot use radiative transfer without gas, --hydro must be "
        "chosen\n");
  }
  if (with_rt && !with_stars) {
    error(
        "Error: Cannot use radiative transfer without stars, --stars must be "
        "chosen\n");
  }
  if (with_rt && !with_feedback) {
    error(
        "Error: Cannot use radiative transfer without --feedback "
        "(even if configured --with-feedback=none)");
  }

#ifdef WITH_MPI
  /* Temporary, this will be removed in due time */
  error("Error: Cannot use radiative transfer with MPI\n");
#endif
#endif /* idfef RT_NONE */

#ifdef SINK_NONE
  if (with_sink) {
    error("Running with sink particles but compiled without them!");
  }
#endif

/* Let's pin the main thread, now we know if affinity will be used. */
#if defined(HAVE_SETAFFINITY) && defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
  if (with_aff &&
      ((ENGINE_POLICY)&engine_policy_setaffinity) == engine_policy_setaffinity)
    engine_pin();
#endif

  /* Genesis 1.1: And then, there was time ! */
  clocks_set_cpufreq(cpufreq);

  /* Are we running with gravity */
  const int with_gravity = (with_self_gravity || with_external_gravity);

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

/* Do we have debugging checks ? */
#ifdef SWIFT_USE_NAIVE_INTERACTIONS
  if (myrank == 0)
    message(
        "WARNING: Naive cell interactions activated. Code will be slower !");
#endif

/* Do we have gravity accuracy checks ? */
#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  if (myrank == 0)
    message(
        "WARNING: Checking 1/%d of all gpart for gravity accuracy. Code will "
        "be slower !",
        SWIFT_GRAVITY_FORCE_CHECKS);
#endif

/* Do we have hydro accuracy checks ? */
#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  if (myrank == 0)
    message(
        "WARNING: Checking 1/%d of all part for hydro accuracy. Code will be "
        "slower !",
        SWIFT_HYDRO_DENSITY_CHECKS);
#endif

  /* Do we choke on FP-exceptions ? */
  if (with_fp_exceptions) {
#ifdef HAVE_FE_ENABLE_EXCEPT
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
    if (myrank == 0)
      message("WARNING: Floating point exceptions will be reported.");
  }

/* Do we have slow barriers? */
#ifndef HAVE_PTHREAD_BARRIERS
  if (myrank == 0)
    message("WARNING: Non-optimal thread barriers are being used.");
#endif

  /* How large are the parts? */
  if (myrank == 0) {
    message("sizeof(part)        is %4zi bytes.", sizeof(struct part));
    message("sizeof(xpart)       is %4zi bytes.", sizeof(struct xpart));
    message("sizeof(sink)        is %4zi bytes.", sizeof(struct sink));
    message("sizeof(spart)       is %4zi bytes.", sizeof(struct spart));
    message("sizeof(bpart)       is %4zi bytes.", sizeof(struct bpart));
    message("sizeof(gpart)       is %4zi bytes.", sizeof(struct gpart));
    message("sizeof(multipole)   is %4zi bytes.", sizeof(struct multipole));
    message("sizeof(grav_tensor) is %4zi bytes.", sizeof(struct grav_tensor));
    message("sizeof(task)        is %4zi bytes.", sizeof(struct task));
    message("sizeof(cell)        is %4zi bytes.", sizeof(struct cell));
  }

  /* Read the parameter file */
  struct swift_params *params =
      (struct swift_params *)malloc(sizeof(struct swift_params));
  if (params == NULL) error("Error allocating memory for the parameter file.");
  if (myrank == 0) {
    message("Reading runtime parameters from file '%s'", param_filename);
    parser_read_file(param_filename, params);

    /* Handle any command-line overrides. */
    if (cmdps.nparam > 0) {
      message(
          "Overwriting values read from the YAML file with command-line "
          "values.");
      for (int k = 0; k < cmdps.nparam; k++)
        parser_set_param(params, cmdps.param[k]);
    }
  }

#ifdef WITH_MPI
  /* Broadcast the parameter file */
  MPI_Bcast(params, sizeof(struct swift_params), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  /* Read the provided output selection file, if available. Best to
   * do this after broadcasting the parameters as there may be code in this
   * function that is repeated on each node based on the parameter file. */

  struct output_options *output_options =
      (struct output_options *)malloc(sizeof(struct output_options));
  output_options_init(params, myrank, output_options);

  /* Temporary early aborts for modes not supported over MPI. */
#ifdef WITH_MPI
  if (with_mpole_reconstruction && nr_nodes > 1)
    error("Cannot reconstruct m-poles every step over MPI (yet).");
#endif

    /* Temporary early aborts for modes not supported with hand-vec. */
#if defined(WITH_VECTORIZATION) && defined(GADGET2_SPH) && \
    !defined(CHEMISTRY_NONE)
  error(
      "Cannot run with chemistry and hand-vectorization (yet). "
      "Use --disable-hand-vec at configure time.");
#endif

  /* Check that we can write the snapshots by testing if the output
   * directory exists and is searchable and writable. */
  char basename[PARSER_MAX_LINE_SIZE];
  parser_get_param_string(params, "Snapshots:basename", basename);
  const char *dirp = dirname(basename);
  if (access(dirp, W_OK | X_OK) != 0) {
    error("Cannot write snapshots in directory %s (%s)", dirp, strerror(errno));
  }

  /* Check that we can write the structure finding catalogues by testing if the
   * output directory exists and is searchable and writable. */
  if (with_structure_finding) {
    char stfbasename[PARSER_MAX_LINE_SIZE];
    parser_get_param_string(params, "StructureFinding:basename", stfbasename);
    const char *stfdirp = dirname(stfbasename);
    if (access(stfdirp, W_OK | X_OK) != 0) {
      error("Cannot write stf catalogues in directory %s (%s)", stfdirp,
            strerror(errno));
    }
  }

  /* Check that we can write the line of sight files by testing if the
   * output directory exists and is searchable and writable. */
  if (with_line_of_sight) {
    char losbasename[PARSER_MAX_LINE_SIZE];
    parser_get_param_string(params, "LineOfSight:basename", losbasename);
    const char *losdirp = dirname(losbasename);
    if (access(losdirp, W_OK | X_OK) != 0) {
      error("Cannot write line of sight in directory %s (%s)", losdirp,
            strerror(errno));
    }
  }

  /* Prepare the domain decomposition scheme */
  struct repartition reparttype;
#ifdef WITH_MPI
  struct partition initial_partition;
  partition_init(&initial_partition, &reparttype, params, nr_nodes);

  /* Let's report what we did */
  if (myrank == 0) {
#if defined(HAVE_PARMETIS)
    if (reparttype.usemetis)
      message("Using METIS serial partitioning:");
    else
      message("Using ParMETIS partitioning:");
#elif defined(HAVE_METIS)
    message("Using METIS serial partitioning:");
#else
    message("Non-METIS partitioning:");
#endif
    message("  initial partitioning: %s",
            initial_partition_name[initial_partition.type]);
    if (initial_partition.type == INITPART_GRID)
      message("    grid set to [ %i %i %i ].", initial_partition.grid[0],
              initial_partition.grid[1], initial_partition.grid[2]);
    message("  repartitioning: %s", repartition_name[reparttype.type]);
  }
#endif

  /* Common variables for restart and IC sections. */
  int clean_smoothing_length_values = 0;
  int flag_entropy_ICs = 0;

  /* Work out where we will read and write restart files. */
  char restart_dir[PARSER_MAX_LINE_SIZE];
  parser_get_opt_param_string(params, "Restarts:subdir", restart_dir,
                              "restart");

  /* The directory must exist. */
  if (myrank == 0) {
    if (access(restart_dir, W_OK | X_OK) != 0) {
      if (restart) {
        error("Cannot restart as no restart subdirectory: %s (%s)", restart_dir,
              strerror(errno));
      } else {
        if (mkdir(restart_dir, 0777) != 0)
          error("Failed to create restart directory: %s (%s)", restart_dir,
                strerror(errno));
      }
    }
  }

  /* Basename for any restart files. */
  char restart_name[PARSER_MAX_LINE_SIZE];
  parser_get_opt_param_string(params, "Restarts:basename", restart_name,
                              "swift");

  /* How often to check for the stop file and dump restarts and exit the
   * application. */
  const int restart_stop_steps =
      parser_get_opt_param_int(params, "Restarts:stop_steps", 100);

  /* Get the maximal wall-clock time of this run */
  const float restart_max_hours_runtime =
      parser_get_opt_param_float(params, "Restarts:max_run_time", FLT_MAX);

  /* Do we want to resubmit when we hit the limit? */
  const int resubmit_after_max_hours =
      parser_get_opt_param_int(params, "Restarts:resubmit_on_exit", 0);

  /* What command should we run to resubmit at the end? */
  char resubmit_command[PARSER_MAX_LINE_SIZE];
  if (resubmit_after_max_hours)
    parser_get_param_string(params, "Restarts:resubmit_command",
                            resubmit_command);

  /* If restarting, look for the restart files. */
  if (restart) {

    /* Attempting a restart. */
    char **restart_files = NULL;
    int restart_nfiles = 0;

    if (myrank == 0) {
      message("Restarting SWIFT");

      /* Locate the restart files. */
      restart_files =
          restart_locate(restart_dir, restart_name, &restart_nfiles);
      if (restart_nfiles == 0)
        error("Failed to locate any restart files in %s", restart_dir);

      /* We need one file per rank. */
      if (restart_nfiles != nr_nodes)
        error("Incorrect number of restart files, expected %d found %d",
              nr_nodes, restart_nfiles);

      if (verbose > 0)
        for (int i = 0; i < restart_nfiles; i++)
          message("found restart file: %s", restart_files[i]);
    }

#ifdef WITH_MPI
    /* Distribute the restart files, need one for each rank. */
    if (myrank == 0) {

      for (int i = 1; i < nr_nodes; i++) {
        strcpy(restart_file, restart_files[i]);
        MPI_Send(restart_file, 200, MPI_BYTE, i, 0, MPI_COMM_WORLD);
      }

      /* Keep local file. */
      strcpy(restart_file, restart_files[0]);

      /* Finished with the list. */
      restart_locate_free(restart_nfiles, restart_files);

    } else {
      MPI_Recv(restart_file, 200, MPI_BYTE, 0, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
    }
    if (verbose > 1) message("local restart file = %s", restart_file);
#else

    /* Just one restart file. */
    strcpy(restart_file, restart_files[0]);

    /* Finished with the list. */
    restart_locate_free(1, restart_files);
#endif

    /* Now read it. */
    restart_read(&e, restart_file);

#ifdef WITH_MPI
    integertime_t min_ti_current = e.ti_current;
    integertime_t max_ti_current = e.ti_current;

    /* Verify that everyone agrees on the current time */
    MPI_Allreduce(&e.ti_current, &min_ti_current, 1, MPI_LONG_LONG_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&e.ti_current, &max_ti_current, 1, MPI_LONG_LONG_INT, MPI_MAX,
                  MPI_COMM_WORLD);

    if (min_ti_current != max_ti_current) {
      if (myrank == 0)
        message("The restart files don't all contain the same ti_current!");

      for (int i = 0; i < myrank; ++i) {
        if (myrank == i)
          message("MPI rank %d reading file '%s' found an integer time= %lld",
                  myrank, restart_file, e.ti_current);
        MPI_Barrier(MPI_COMM_WORLD);
      }

      if (myrank == 0) error("Aborting");
    }
#endif

    /* And initialize the engine with the space and policies. */
    engine_config(/*restart=*/1, /*fof=*/0, &e, params, nr_nodes, myrank,
                  nr_threads, nr_pool_threads, with_aff, talking, restart_file);

    /* Check if we are already done when given steps on the command-line. */
    if (e.step >= nsteps && nsteps > 0)
      error("Not restarting, already completed %d steps", e.step);

  } else {

    /* Prepare and verify the selection of outputs */
    io_prepare_output_fields(output_options, with_cosmology, with_fof,
                             with_structure_finding, e.verbose);

#if defined(SWIFT_DEBUG_TASKS)
    task_create_name_files("task_labels");
#endif

    /* Not restarting so look for the ICs. */
    /* Initialize unit system and constants */
    units_init_from_params(&us, params, "InternalUnitSystem");
    phys_const_init(&us, params, &prog_const);
    if (myrank == 0) {
      message("Internal unit system: U_M = %e g.", us.UnitMass_in_cgs);
      message("Internal unit system: U_L = %e cm.", us.UnitLength_in_cgs);
      message("Internal unit system: U_t = %e s.", us.UnitTime_in_cgs);
      message("Internal unit system: U_I = %e A.", us.UnitCurrent_in_cgs);
      message("Internal unit system: U_T = %e K.", us.UnitTemperature_in_cgs);
      phys_const_print(&prog_const);
    }

    /* Read particles and space information from ICs */
    char ICfileName[200] = "";
    parser_get_param_string(params, "InitialConditions:file_name", ICfileName);
    const int periodic =
        parser_get_param_int(params, "InitialConditions:periodic");
    const int replicate =
        parser_get_opt_param_int(params, "InitialConditions:replicate", 1);
    clean_smoothing_length_values = parser_get_opt_param_int(
        params, "InitialConditions:cleanup_smoothing_lengths", 0);
    const int cleanup_h = parser_get_opt_param_int(
        params, "InitialConditions:cleanup_h_factors", 0);
    const int cleanup_sqrt_a = parser_get_opt_param_int(
        params, "InitialConditions:cleanup_velocity_factors", 0);
    const int generate_gas_in_ics = parser_get_opt_param_int(
        params, "InitialConditions:generate_gas_in_ics", 0);
    const int remap_ids =
        parser_get_opt_param_int(params, "InitialConditions:remap_ids", 0);

    /* Initialise the cosmology */
    if (with_cosmology)
      cosmology_init(params, &us, &prog_const, &cosmo);
    else
      cosmology_init_no_cosmo(&cosmo);
    if (myrank == 0 && with_cosmology) cosmology_print(&cosmo);

    if (with_hydro) {
#ifdef NONE_SPH
      error("Can't run with hydro when compiled without a hydro model!");
#endif
    }
    if (with_stars) {
#ifdef STARS_NONE
      error("Can't run with stars when compiled without a stellar model!");
#endif
    }
    if (with_black_holes) {
#ifdef BLACK_HOLES_NONE
      error(
          "Can't run with black holes when compiled without a black hole "
          "model!");
#endif
    }

    /* Initialise the hydro properties */
    if (with_hydro)
      hydro_props_init(&hydro_properties, &prog_const, &us, params);
    else
      bzero(&hydro_properties, sizeof(struct hydro_props));

    /* Initialise the equation of state */
    if (with_hydro)
      eos_init(&eos, &prog_const, &us, params);
    else
      bzero(&eos, sizeof(struct eos_parameters));

    /* Initialise the entropy floor */
    if (with_hydro)
      entropy_floor_init(&entropy_floor, &prog_const, &us, &hydro_properties,
                         params);
    else
      bzero(&entropy_floor, sizeof(struct entropy_floor_properties));

    /* Initialise the pressure floor */
    if (with_hydro)
      pressure_floor_init(&pressure_floor_props, &prog_const, &us,
                          &hydro_properties, params);
    else
      bzero(&pressure_floor_props, sizeof(struct pressure_floor_properties));

    /* Initialise the stars properties */
    if (with_stars)
      stars_props_init(&stars_properties, &prog_const, &us, params,
                       &hydro_properties, &cosmo);
    else
      bzero(&stars_properties, sizeof(struct stars_props));

    /* Initialise the feedback properties */
    if (with_feedback) {
#ifdef FEEDBACK_NONE
      if (!with_rt) /* allow swift to run without feedback model if RT is on. */
        error("ERROR: Running with feedback but compiled without it.");
#endif
      feedback_props_init(&feedback_properties, &prog_const, &us, params,
                          &hydro_properties, &cosmo);
    } else
      bzero(&feedback_properties, sizeof(struct feedback_props));

    /* Initialize RT properties */
    if (with_rt) {
      if (hydro_properties.particle_splitting)
        error("Can't run with RT and particle splitting as of yet.");
      rt_props_init(&rt_properties, &prog_const, &us, params, &cosmo);
    } else
      bzero(&rt_properties, sizeof(struct rt_props));

    /* Initialise the black holes properties */
    if (with_black_holes) {
#ifdef BLACK_HOLES_NONE
      error("ERROR: Running with black_holes but compiled without it.");
#endif
      black_holes_props_init(&black_holes_properties, &prog_const, &us, params,
                             &hydro_properties, &cosmo);
    } else
      bzero(&black_holes_properties, sizeof(struct black_holes_props));

    /* Initialise the sink properties */
    if (with_sink) {
      sink_props_init(&sink_properties, &prog_const, &us, params, &cosmo);
    } else
      bzero(&sink_properties, sizeof(struct sink_props));

      /* Initialise the cooling function properties */
#ifdef COOLING_NONE
    if (with_cooling) {
      error(
          "ERROR: Running with cooling calculation"
          " but compiled without it.");
    }
#endif
    bzero(&cooling_func, sizeof(struct cooling_function_data));
    if (with_cooling || with_temperature) {
      cooling_init(params, &us, &prog_const, &hydro_properties, &cooling_func);
    }
    if (myrank == 0) cooling_print(&cooling_func);

    /* Initialise the star formation law and its properties */
    bzero(&starform, sizeof(struct star_formation));
    if (with_star_formation) {
#ifdef STAR_FORMATION_NONE
      error("ERROR: Running with star formation but compiled without it!");
#endif
      starformation_init(params, &prog_const, &us, &hydro_properties, &cosmo,
                         &entropy_floor, &starform);
    }
    if (with_star_formation && myrank == 0) starformation_print(&starform);

    /* Initialise the chemistry */
    bzero(&chemistry, sizeof(struct chemistry_global_data));
    chemistry_init(params, &us, &prog_const, &chemistry);
    if (myrank == 0) chemistry_print(&chemistry);

    /* Initialise the FOF properties */
    bzero(&fof_properties, sizeof(struct fof_props));
#ifdef WITH_FOF
    if (with_fof) {
      fof_init(&fof_properties, params, &prog_const, &us, /*stand-alone=*/0);
      if (fof_properties.seed_black_holes_enabled && !with_black_holes) {
        if (myrank == 0)
          printf(
              "Error: Cannot perform FOF seeding without black holes being in "
              "use\n");
        return 1;
      }
    }
#endif

    /* Be verbose about what happens next */
    if (myrank == 0) message("Reading ICs from file '%s'", ICfileName);
    if (myrank == 0 && cleanup_h)
      message("Cleaning up h-factors (h=%f)", cosmo.h);
    if (myrank == 0 && cleanup_sqrt_a)
      message("Cleaning up a-factors from velocity (a=%f)", cosmo.a);
    fflush(stdout);

    /* Get ready to read particles of all kinds */
    size_t Ngas = 0, Ngpart = 0, Ngpart_background = 0, Nnupart = 0;
    size_t Nspart = 0, Nbpart = 0, Nsink = 0;
    double dim[3] = {0., 0., 0.};

    /* Prepare struct to store metadata from ICs */
    ic_info_init(&ics_metadata, params);

    if (myrank == 0) clocks_gettime(&tic);
#if defined(HAVE_HDF5)
#if defined(WITH_MPI)
#if defined(HAVE_PARALLEL_HDF5)
    read_ic_parallel(ICfileName, &us, dim, &parts, &gparts, &sinks, &sparts,
                     &bparts, &Ngas, &Ngpart, &Ngpart_background, &Nnupart,
                     &Nsink, &Nspart, &Nbpart, &flag_entropy_ICs, with_hydro,
                     with_gravity, with_sink, with_stars, with_black_holes,
                     with_cosmology, cleanup_h, cleanup_sqrt_a, cosmo.h,
                     cosmo.a, myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL,
                     nr_threads, dry_run, remap_ids, &ics_metadata);
#else
    read_ic_serial(ICfileName, &us, dim, &parts, &gparts, &sinks, &sparts,
                   &bparts, &Ngas, &Ngpart, &Ngpart_background, &Nnupart,
                   &Nsink, &Nspart, &Nbpart, &flag_entropy_ICs, with_hydro,
                   with_gravity, with_sink, with_stars, with_black_holes,
                   with_cosmology, cleanup_h, cleanup_sqrt_a, cosmo.h, cosmo.a,
                   myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL, nr_threads,
                   dry_run, remap_ids, &ics_metadata);
#endif
#else
    read_ic_single(ICfileName, &us, dim, &parts, &gparts, &sinks, &sparts,
                   &bparts, &Ngas, &Ngpart, &Ngpart_background, &Nnupart,
                   &Nsink, &Nspart, &Nbpart, &flag_entropy_ICs, with_hydro,
                   with_gravity, with_sink, with_stars, with_black_holes,
                   with_cosmology, cleanup_h, cleanup_sqrt_a, cosmo.h, cosmo.a,
                   nr_threads, dry_run, remap_ids, &ics_metadata);
#endif
#endif

    if (myrank == 0) {
      clocks_gettime(&toc);
      message("Reading initial conditions took %.3f %s.",
              clocks_diff(&tic, &toc), clocks_getunit());
      fflush(stdout);
    }

    /* Some checks that we are not doing something stupid */
    if (generate_gas_in_ics && flag_entropy_ICs)
      error("Can't generate gas if the entropy flag is set in the ICs.");
    if (generate_gas_in_ics && !with_cosmology)
      error("Can't generate gas if the run is not cosmological.");

#ifdef SWIFT_DEBUG_CHECKS
    /* Check once and for all that we don't have unwanted links */
    for (size_t k = 0; k < Ngpart; ++k)
      if (!dry_run && gparts[k].id_or_neg_offset == 0 &&
          (gparts[k].type == swift_type_dark_matter ||
           gparts[k].type == swift_type_dark_matter_background))
        error("SWIFT does not allow the ID 0.");
    if (!with_stars && !dry_run) {
      for (size_t k = 0; k < Ngpart; ++k)
        if (gparts[k].type == swift_type_stars) error("Linking problem");
    }
    if (!with_black_holes && !dry_run) {
      for (size_t k = 0; k < Ngpart; ++k)
        if (gparts[k].type == swift_type_black_hole) error("Linking problem");
    }
    if (!with_hydro && !dry_run) {
      for (size_t k = 0; k < Ngpart; ++k)
        if (gparts[k].type == swift_type_gas) error("Linking problem");
    }
    if (!with_sink && !dry_run) {
      for (size_t k = 0; k < Ngpart; ++k)
        if (gparts[k].type == swift_type_sink) error("Linking problem");
    }

    /* Check that the other links are correctly set */
    if (!dry_run)
      part_verify_links(parts, gparts, sinks, sparts, bparts, Ngas, Ngpart,
                        Nsink, Nspart, Nbpart, /*verbose=*/1);
#endif

    /* Get the total number of particles across all nodes. */
    long long N_total[swift_type_count + 1] = {0};
    long long Nbaryons = Ngas + Nspart + Nbpart + Nsink;
#if defined(WITH_MPI)
    long long N_long[swift_type_count + 1] = {0};
    N_long[swift_type_gas] = Ngas;
    N_long[swift_type_dark_matter_background] = Ngpart_background;
    N_long[swift_type_sink] = Nsink;
    N_long[swift_type_stars] = Nspart;
    N_long[swift_type_black_hole] = Nbpart;
    N_long[swift_type_neutrino] = Nnupart;
    N_long[swift_type_count] = Ngpart;
    N_long[swift_type_dark_matter] =
        with_gravity ? Ngpart - Ngpart_background - Nbaryons - Nnupart : 0;

    MPI_Allreduce(&N_long, &N_total, swift_type_count + 1, MPI_LONG_LONG_INT,
                  MPI_SUM, MPI_COMM_WORLD);
#else
    N_total[swift_type_gas] = Ngas;
    N_total[swift_type_dark_matter_background] = Ngpart_background;
    N_total[swift_type_sink] = Nsink;
    N_total[swift_type_stars] = Nspart;
    N_total[swift_type_black_hole] = Nbpart;
    N_total[swift_type_neutrino] = Nnupart;
    N_total[swift_type_count] = Ngpart;
    N_total[swift_type_dark_matter] =
        with_gravity ? Ngpart - Ngpart_background - Nbaryons - Nnupart : 0;
#endif

    if (myrank == 0)
      message(
          "Read %lld gas particles, %lld sink particles, %lld star particles, "
          "%lld black hole particles, %lld DM particles, %lld DM background "
          "particles, and %lld neutrino DM particles from the ICs.",
          N_total[swift_type_gas], N_total[swift_type_sink],
          N_total[swift_type_stars], N_total[swift_type_black_hole],
          N_total[swift_type_dark_matter],
          N_total[swift_type_dark_matter_background],
          N_total[swift_type_neutrino]);

    const int with_DM_particles = N_total[swift_type_dark_matter] > 0;
    const int with_baryon_particles =
        (N_total[swift_type_gas] + N_total[swift_type_stars] +
             N_total[swift_type_black_hole] + N_total[swift_type_sink] >
         0) ||
        (with_DM_particles && generate_gas_in_ics);

    /* Do we have background DM particles? */
    const int with_DM_background_particles =
        N_total[swift_type_dark_matter_background] > 0;

    /* Do we have neutrino DM particles? */
    const int with_neutrinos = N_total[swift_type_neutrino] > 0;

    /* Initialise the neutrino properties if we have neutrino particles */
    bzero(&neutrino_properties, sizeof(struct neutrino_props));
    if (with_neutrinos)
      neutrino_props_init(&neutrino_properties, &prog_const, &us, params,
                          &cosmo);

    /* Initialize the space with these data. */
    if (myrank == 0) clocks_gettime(&tic);
    space_init(&s, params, &cosmo, dim, &hydro_properties, parts, gparts, sinks,
               sparts, bparts, Ngas, Ngpart, Nsink, Nspart, Nbpart, Nnupart,
               periodic, replicate, remap_ids, generate_gas_in_ics, with_hydro,
               with_self_gravity, with_star_formation, with_sink,
               with_DM_background_particles, with_neutrinos, talking, dry_run,
               nr_nodes);

    /* Initialise the line of sight properties. */
    if (with_line_of_sight) los_init(s.dim, &los_properties, params);

    if (myrank == 0) {
      clocks_gettime(&toc);
      message("space_init took %.3f %s.", clocks_diff(&tic, &toc),
              clocks_getunit());
      fflush(stdout);
    }

    /* Initialise the gravity properties */
    bzero(&gravity_properties, sizeof(struct gravity_props));
    if (with_self_gravity)
      gravity_props_init(
          &gravity_properties, params, &prog_const, &cosmo, with_cosmology,
          with_external_gravity, with_baryon_particles, with_DM_particles,
          with_neutrinos, with_DM_background_particles, periodic, s.dim);

    /* Initialise the external potential properties */
    bzero(&potential, sizeof(struct external_potential));
    if (with_external_gravity)
      potential_init(params, &prog_const, &us, &s, &potential);
    if (myrank == 0) potential_print(&potential);

    /* Initialise the long-range gravity mesh */
    if (with_self_gravity && periodic) {
#ifdef HAVE_FFTW
      pm_mesh_init(&mesh, &gravity_properties, s.dim, nr_threads);
#else
      /* Need the FFTW library if periodic and self gravity. */
      error(
          "No FFTW library found. Cannot compute periodic long-range forces.");
#endif
    } else {
      pm_mesh_init_no_mesh(&mesh, s.dim);
    }

    /* Also update the total counts (in case of changes due to replication) */
    Nbaryons = s.nr_parts + s.nr_sparts + s.nr_bparts + s.nr_sinks;
    Nnupart = s.nr_nuparts;
#if defined(WITH_MPI)
    N_long[swift_type_gas] = s.nr_parts;
    N_long[swift_type_dark_matter] =
        with_gravity ? s.nr_gparts - Ngpart_background - Nbaryons - Nnupart : 0;
    N_long[swift_type_count] = s.nr_gparts;
    N_long[swift_type_stars] = s.nr_sparts;
    N_long[swift_type_sink] = s.nr_sinks;
    N_long[swift_type_black_hole] = s.nr_bparts;
    N_long[swift_type_neutrino] = s.nr_nuparts;
    MPI_Allreduce(&N_long, &N_total, swift_type_count + 1, MPI_LONG_LONG_INT,
                  MPI_SUM, MPI_COMM_WORLD);
#else
    N_total[swift_type_gas] = s.nr_parts;
    N_total[swift_type_dark_matter] =
        with_gravity ? s.nr_gparts - Ngpart_background - Nbaryons - Nnupart : 0;
    N_total[swift_type_count] = s.nr_gparts;
    N_total[swift_type_stars] = s.nr_sparts;
    N_total[swift_type_sink] = s.nr_sinks;
    N_total[swift_type_black_hole] = s.nr_bparts;
    N_total[swift_type_neutrino] = s.nr_nuparts;
#endif

    /* Say a few nice things about the space we just created. */
    if (myrank == 0) {
      message("space dimensions are [ %.3f %.3f %.3f ].", s.dim[0], s.dim[1],
              s.dim[2]);
      message("space %s periodic.", s.periodic ? "is" : "isn't");
      message("highest-level cell dimensions are [ %i %i %i ].", s.cdim[0],
              s.cdim[1], s.cdim[2]);
      message("%zi parts in %i cells.", s.nr_parts, s.tot_cells);
      message("%zi gparts in %i cells.", s.nr_gparts, s.tot_cells);
      message("%zi sinks in %i cells.", s.nr_sinks, s.tot_cells);
      message("%zi sparts in %i cells.", s.nr_sparts, s.tot_cells);
      message("%zi bparts in %i cells.", s.nr_bparts, s.tot_cells);
      message("maximum depth is %d.", s.maxdepth);
      fflush(stdout);
    }

    /* Verify that we are not using basic modes incorrectly */
    if (with_hydro && N_total[swift_type_gas] == 0) {
      error(
          "ERROR: Running with hydrodynamics but no gas particles found in the "
          "ICs!");
    }
    if (with_gravity && N_total[swift_type_count] == 0) {
      error(
          "ERROR: Running with gravity but no gravity particles found in "
          "the ICs!");
    }

    /* Verify that each particle is in its proper cell. */
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

    /* Construct the engine policy */
    int engine_policies = ENGINE_POLICY | engine_policy_steal;
    if (with_drift_all) engine_policies |= engine_policy_drift_all;
    if (with_mpole_reconstruction)
      engine_policies |= engine_policy_reconstruct_mpoles;
    if (with_hydro) engine_policies |= engine_policy_hydro;
    if (with_self_gravity) engine_policies |= engine_policy_self_gravity;
    if (with_external_gravity)
      engine_policies |= engine_policy_external_gravity;
    if (with_cosmology) engine_policies |= engine_policy_cosmology;
    if (with_temperature) engine_policies |= engine_policy_temperature;
    if (with_timestep_limiter)
      engine_policies |= engine_policy_timestep_limiter;
    if (with_timestep_sync) engine_policies |= engine_policy_timestep_sync;
    if (with_cooling) engine_policies |= engine_policy_cooling;
    if (with_stars) engine_policies |= engine_policy_stars;
    if (with_star_formation) engine_policies |= engine_policy_star_formation;
    if (with_feedback) engine_policies |= engine_policy_feedback;
    if (with_black_holes) engine_policies |= engine_policy_black_holes;
    if (with_structure_finding)
      engine_policies |= engine_policy_structure_finding;
    if (with_fof) engine_policies |= engine_policy_fof;
    if (with_csds) engine_policies |= engine_policy_csds;
    if (with_line_of_sight) engine_policies |= engine_policy_line_of_sight;
    if (with_sink) engine_policies |= engine_policy_sinks;
    if (with_rt) engine_policies |= engine_policy_rt;

    /* Initialize the engine with the space and policies. */
    engine_init(&e, &s, params, output_options, N_total[swift_type_gas],
                N_total[swift_type_count], N_total[swift_type_sink],
                N_total[swift_type_stars], N_total[swift_type_black_hole],
                N_total[swift_type_dark_matter_background],
                N_total[swift_type_neutrino], engine_policies, talking,
                &reparttype, &us, &prog_const, &cosmo, &hydro_properties,
                &entropy_floor, &gravity_properties, &stars_properties,
                &black_holes_properties, &sink_properties, &neutrino_properties,
                &feedback_properties, &rt_properties, &mesh, &potential,
                &cooling_func, &starform, &chemistry, &fof_properties,
                &los_properties, &ics_metadata);
    engine_config(/*restart=*/0, /*fof=*/0, &e, params, nr_nodes, myrank,
                  nr_threads, nr_pool_threads, with_aff, talking, restart_file);

    /* Compute some stats for the star formation */
    if (with_star_formation) {
      star_formation_first_init_stats(&starform, &e);
    }

    /* Get some info to the user. */
    if (myrank == 0) {
      const long long N_DM = N_total[swift_type_dark_matter] +
                             N_total[swift_type_dark_matter_background];
      message(
          "Running on %lld gas particles, %lld sink particles, %lld stars "
          "particles %lld black hole particles, %lld neutrino particles, and "
          "%lld DM particles (%lld gravity particles)",
          N_total[swift_type_gas], N_total[swift_type_sink],
          N_total[swift_type_stars], N_total[swift_type_black_hole],
          N_total[swift_type_neutrino], N_DM, N_total[swift_type_count]);
      message(
          "from t=%.3e until t=%.3e with %d ranks, %d threads / rank and %d "
          "task queues / rank (dt_min=%.3e, dt_max=%.3e)...",
          e.time_begin, e.time_end, nr_nodes, e.nr_threads, e.sched.nr_queues,
          e.dt_min, e.dt_max);
      fflush(stdout);
    }
  }

  /* Time to say good-bye if this was not a serious run. */
  if (dry_run) {
#ifdef WITH_MPI
    if ((res = MPI_Finalize()) != MPI_SUCCESS)
      error("call to MPI_Finalize failed with error %i.", res);
#endif
    if (myrank == 0)
      message("Time integration ready to start. End of dry-run.");
    engine_clean(&e, /*fof=*/0, /*restart=*/0);
    free(params);
    return 0;
  }

/* Initialise the table of Ewald corrections for the gravity checks */
#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  if (s.periodic) gravity_exact_force_ewald_init(e.s->dim[0]);
#endif

  if (!restart) {

#ifdef WITH_MPI
    /* Split the space. */
    engine_split(&e, &initial_partition);
#endif

    /* Initialise the particles */
    engine_init_particles(&e, flag_entropy_ICs, clean_smoothing_length_values);

    /* Check that the matter content matches the cosmology given in the
     * parameter file. */
    if (with_cosmology && with_self_gravity && !dry_run) {
      const int check_neutrinos = !neutrino_properties.use_delta_f;
      space_check_cosmology(&s, &cosmo, with_hydro, myrank, check_neutrinos);
      neutrino_check_cosmology(&s, &cosmo, &prog_const, params,
                               &neutrino_properties, myrank, verbose);
    }

    /* Write the state of the system before starting time integration. */
#ifdef WITH_CSDS
    if (e.policy & engine_policy_csds) {
      csds_log_all_particles(e.csds, &e, csds_flag_create);
    }
#endif
    /* Dump initial state snapshot, if not working with an output list */
    if (!e.output_list_snapshots) {

      /* Run FoF first, if we're adding FoF info to the snapshot */
      if (with_fof && e.snapshot_invoke_fof) {
        engine_fof(&e, /*dump_results=*/1, /*dump_debug=*/0,
                   /*seed_black_holes=*/0, /*buffers allocated=*/1);
      }

      engine_dump_snapshot(&e);
    }

    /* Dump initial state statistics, if not working with an output list */
    if (!e.output_list_stats) engine_print_stats(&e);

    /* Is there a dump before the end of the first time-step? */
    engine_check_for_dumps(&e);
  }

  /* Legend */
  if (myrank == 0) {
    printf(
        "# %6s %14s %12s %12s %14s %9s %12s %12s %12s %12s %12s %16s [%s] "
        "%6s %s [%s] \n",
        "Step", "Time", "Scale-factor", "Redshift", "Time-step", "Time-bins",
        "Updates", "g-Updates", "s-Updates", "sink-Updates", "b-Updates",
        "Wall-clock time", clocks_getunit(), "Props", "Dead time",
        clocks_getunit());
    fflush(stdout);
  }

  /* File for the timers */
  if (with_verbose_timers) timers_open_file(myrank);

  /* Create a name for restart file of this rank. */
  if (restart_genname(restart_dir, restart_name, e.nodeID, restart_file, 200) !=
      0)
    error("Failed to generate restart filename");

  /* dump the parameters as used. */
  if (!restart && myrank == 0) {

    /* used parameters */
    parser_write_params_to_file(params, "used_parameters.yml", /*used=*/1);
    /* unused parameters */
    parser_write_params_to_file(params, "unused_parameters.yml", /*used=*/0);
  }

  /* Dump memory use report if collected for the 0 step. */
#ifdef SWIFT_MEMUSE_REPORTS
  {
    char dumpfile[40];
#ifdef WITH_MPI
    snprintf(dumpfile, 40, "memuse_report-rank%d-step%d.dat", engine_rank, 0);
#else
    snprintf(dumpfile, 40, "memuse_report-step%d.dat", 0);
#endif  // WITH_MPI
    memuse_log_dump(dumpfile);
  }
#endif

  /* Dump MPI requests if collected. */
#if defined(SWIFT_MPIUSE_REPORTS) && defined(WITH_MPI)
  {
    char dumpfile[40];
    snprintf(dumpfile, 40, "mpiuse_report-rank%d-step%d.dat", engine_rank, 0);
    mpiuse_log_dump(dumpfile, clocks_start_ticks);
  }
#endif

  /* Main simulation loop */
  /* ==================== */
  int force_stop = 0, resubmit = 0;
  for (int j = 0; !engine_is_done(&e) && e.step - 1 != nsteps && !force_stop;
       j++) {

    /* Reset timers */
    timers_reset_all();

    /* Take a step. */
    engine_step(&e);

    /* Print the timers. */
    if (with_verbose_timers) timers_print(e.step);

    /* Every so often allow the user to stop the application and dump the
     * restart files. */
    if (j % restart_stop_steps == 0) {
      force_stop = restart_stop_now(restart_dir, 0);
      if (myrank == 0 && force_stop)
        message("Forcing application exit, dumping restart files...");
    }

    /* Did we exceed the maximal runtime? */
    if (e.runtime > restart_max_hours_runtime) {
      force_stop = 1;
      message("Runtime limit reached, dumping restart files...");
      if (resubmit_after_max_hours) resubmit = 1;
    }

    /* Also if using nsteps to exit, will not have saved any restarts on exit,
     * make sure we do that (useful in testing only). */
    if (force_stop || (e.restart_onexit && e.step - 1 == nsteps))
      engine_dump_restarts(&e, 0, 1);

    /* Dump the task data using the given frequency. */
    if (dump_tasks && (dump_tasks == 1 || j % dump_tasks == 1)) {
#ifdef SWIFT_DEBUG_TASKS
      if (dump_tasks_threshold == 0.) task_dump_all(&e, j + 1);
#endif

      /* Generate the task statistics. */
      char dumpfile[40];
      snprintf(dumpfile, 40, "thread_stats-step%d.dat", e.step + 1);
      task_dump_stats(dumpfile, &e, dump_tasks_threshold,
                      /* header = */ 0, /* allranks = */ 1);
    }

#ifdef SWIFT_CELL_GRAPH
    /* Dump the cell data using the given frequency. */
    if (dump_cells && (dump_cells == 1 || j % dump_cells == 1)) {
      space_write_cell_hierarchy(e.s, j + 1);
    }
#endif

    /* Dump memory use report if collected. */
#ifdef SWIFT_MEMUSE_REPORTS
    {
      char dumpfile[40];
#ifdef WITH_MPI
      snprintf(dumpfile, 40, "memuse_report-rank%d-step%d.dat", engine_rank,
               j + 1);
#else
      snprintf(dumpfile, 40, "memuse_report-step%d.dat", e.step + 1);
#endif  // WITH_MPI
      memuse_log_dump(dumpfile);
    }
#endif

    /* Dump MPI requests if collected. */
#if defined(SWIFT_MPIUSE_REPORTS) && defined(WITH_MPI)
    {
      char dumpfile[80];
      snprintf(dumpfile, 80, "mpiuse_report-rank%d-step%d.dat", engine_rank,
               j + 1);
      mpiuse_log_dump(dumpfile, e.tic_step);
    }
#endif  // WITH_MPI

#ifdef SWIFT_DEBUG_THREADPOOL
    /* Dump the task data using the given frequency. */
    if (dump_threadpool && (dump_threadpool == 1 || j % dump_threadpool == 1)) {
      char dumpfile[80];
#ifdef WITH_MPI
      snprintf(dumpfile, 80, "threadpool_info-rank%d-step%d.dat", engine_rank,
               j + 1);
#else
      snprintf(dumpfile, 80, "threadpool_info-step%d.dat", j + 1);
#endif  // WITH_MPI
      threadpool_dump_log(&e.threadpool, dumpfile, 1);
    } else {
      threadpool_reset_log(&e.threadpool);
    }
#endif  // SWIFT_DEBUG_THREADPOOL
  }

  /* Write final time information */
  if (myrank == 0) {

    const double dead_time = e.global_deadtime / (nr_nodes * e.nr_threads);

    /* Print some information to the screen */
    printf(
        "  %6d %14e %12.7f %12.7f %14e %4d %4d %12lld %12lld %12lld %12lld "
        "%12lld"
        " %21.3f %6d %21.3f\n",
        e.step, e.time, e.cosmology->a, e.cosmology->z, e.time_step,
        e.min_active_bin, e.max_active_bin, e.updates, e.g_updates, e.s_updates,
        e.sink_updates, e.b_updates, e.wallclock_time, e.step_props, dead_time);
    fflush(stdout);

    fprintf(e.file_timesteps,
            "  %6d %14e %12.7f %12.7f %14e %4d %4d %12lld %12lld %12lld %12lld"
            " %12lld %21.3f %6d %21.3f\n",
            e.step, e.time, e.cosmology->a, e.cosmology->z, e.time_step,
            e.min_active_bin, e.max_active_bin, e.updates, e.g_updates,
            e.s_updates, e.sink_updates, e.b_updates, e.wallclock_time,
            e.step_props, dead_time);
    fflush(e.file_timesteps);

    /* Print information to the SFH logger */
    if (e.policy & engine_policy_star_formation) {
      star_formation_logger_write_to_log_file(
          e.sfh_logger, e.time, e.cosmology->a, e.cosmology->z, e.sfh, e.step);
    }
  }

  /* Write final output. */
  if (!force_stop && nsteps == -2) {

    /* Move forward in time */
    e.ti_old = e.ti_current;
    e.ti_current = e.ti_end_min;
    e.max_active_bin = get_max_active_bin(e.ti_end_min);
    e.min_active_bin = get_min_active_bin(e.ti_current, e.ti_old);
    e.step += 1;
    engine_current_step = e.step;

    engine_drift_all(&e, /*drift_mpole=*/0);

    /* Write final statistics? */
    if (e.output_list_stats) {
      if (e.output_list_stats->final_step_dump) engine_print_stats(&e);
    } else {
      engine_print_stats(&e);
    }
#ifdef WITH_CSDS
    if (e.policy & engine_policy_csds) {
      csds_log_all_particles(e.csds, &e, csds_flag_delete);

      /* Write a sentinel timestamp */
      if (e.policy & engine_policy_cosmology) {
        csds_log_timestamp(e.csds, e.ti_current, e.cosmology->a,
                           &e.csds->timestamp_offset);
      } else {
        csds_log_timestamp(e.csds, e.ti_current, e.time,
                           &e.csds->timestamp_offset);
      }
    }
#endif

    /* Write final snapshot? */
    if ((e.output_list_snapshots && e.output_list_snapshots->final_step_dump) ||
        !e.output_list_snapshots) {

      if (with_fof && e.snapshot_invoke_fof) {
        engine_fof(&e, /*dump_results=*/1, /*dump_debug=*/0,
                   /*seed_black_holes=*/0, /*buffers allocated=*/1);
      }

#ifdef HAVE_VELOCIRAPTOR
      if (with_structure_finding && e.snapshot_invoke_stf &&
          !e.stf_this_timestep)
        velociraptor_invoke(&e, /*linked_with_snap=*/1);
#endif
      engine_dump_snapshot(&e);
#ifdef HAVE_VELOCIRAPTOR
      if (with_structure_finding && e.snapshot_invoke_stf &&
          e.s->gpart_group_data)
        swift_free("gpart_group_data", e.s->gpart_group_data);
#endif
    }

    /* Write final stf? */
#ifdef HAVE_VELOCIRAPTOR
    if (with_structure_finding && e.output_list_stf) {
      if (e.output_list_stf->final_step_dump && !e.stf_this_timestep)
        velociraptor_invoke(&e, /*linked_with_snap=*/0);
    }
#endif
  }

  /* Remove the stop file if used. Do this anyway, we could have missed the
   * stop file if normal exit happened first. */
  if (myrank == 0) force_stop = restart_stop_now(restart_dir, 1);

  /* Did we want to run a re-submission command just before dying? */
  if (myrank == 0 && resubmit) {
    message("Running the resubmission command:");
    restart_resubmit(resubmit_command);
    fflush(stdout);
    fflush(stderr);
    message("resubmission command completed.");
  }

  /* Clean everything */
  if (with_verbose_timers) timers_close_file();
  if (with_cosmology) cosmology_clean(e.cosmology);
  if (with_self_gravity) pm_mesh_clean(e.mesh);
  if (with_stars) stars_props_clean(e.stars_properties);
  if (with_cooling || with_temperature) cooling_clean(e.cooling_func);
  if (with_feedback) feedback_clean(e.feedback_props);
  if (with_rt) rt_clean(e.rt_props);
  engine_clean(&e, /*fof=*/0, restart);
  free(params);
  free(output_options);

#ifdef WITH_MPI
  if ((res = MPI_Finalize()) != MPI_SUCCESS)
    error("call to MPI_Finalize failed with error %i.", res);
#endif

  /* Say goodbye. */
  if (myrank == 0) message("done. Bye.");

  /* All is calm, all is bright. */
  return 0;
}
