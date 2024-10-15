/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
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

/* Config parameters. */
#include <config.h>

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
static const char *const fof_usage[] = {
    "fof [options] [[--] param-file]",
    "fof [options] param-file",
    "fof_mpi [options] [[--] param-file]",
    "fof_mpi [options] param-file",
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
  struct cosmology cosmo;
  struct pm_mesh mesh;
  struct gpart *gparts = NULL;
  struct gravity_props gravity_properties;
  struct hydro_props hydro_properties;
  struct fof_props fof_properties;
  struct neutrino_props neutrino_properties;
  struct part *parts = NULL;
  struct phys_const prog_const;
  struct space s;
  struct spart *sparts = NULL;
  struct sink *sinks = NULL;
  struct bpart *bparts = NULL;
  struct unit_system us;
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
  if (myrank == 0) greetings(/*fof=*/1);

  int with_aff = 0;
  int dump_tasks = 0;
  int dump_threadpool = 0;
  int with_fp_exceptions = 0;
  int with_cosmology = 0;
  int with_sinks = 0;
  int with_stars = 0;
  int with_black_holes = 0;
  int with_hydro = 0;
  int verbose = 0;
  int nr_threads = 1;
  char *output_parameters_filename = NULL;
  char *cpufreqarg = NULL;
  char *param_filename = NULL;
  unsigned long long cpufreq = 0;
  struct cmdparams cmdps;
  cmdps.nparam = 0;
  cmdps.param[0] = NULL;
  char *buffer = NULL;

  /* Parse the command-line parameters. */
  struct argparse_option options[] = {
      OPT_HELP(),

      OPT_GROUP("  Friends-of-Friends options:\n"),
      OPT_BOOLEAN('c', "cosmology", &with_cosmology,
                  "Run with cosmological information.", NULL, 0, 0),
      OPT_BOOLEAN(0, "hydro", &with_hydro, "Read gas particles from the ICs.",
                  NULL, 0, 0),
      OPT_BOOLEAN(0, "sinks", &with_sinks, "Read sinks from the ICs.", NULL, 0,
                  0),
      OPT_BOOLEAN(0, "stars", &with_stars, "Read stars from the ICs.", NULL, 0,
                  0),
      OPT_BOOLEAN(0, "black-holes", &with_black_holes,
                  "Read black holes from the ICs.", NULL, 0, 0),

      OPT_GROUP("  Control options:\n"),
      OPT_BOOLEAN('a', "pin", &with_aff,
                  "Pin runners using processor affinity.", NULL, 0, 0),
      OPT_BOOLEAN('e', "fpe", &with_fp_exceptions,
                  "Enable floating-point exceptions (debugging mode).", NULL, 0,
                  0),
      OPT_STRING('f', "cpu-frequency", &cpufreqarg,
                 "Overwrite the CPU "
                 "frequency (Hz) to be used for time measurements.",
                 NULL, 0, 0),
      OPT_STRING('P', "param", &buffer,
                 "Set parameter value, overiding the value read from the "
                 "parameter file. Can be used more than once {sec:par:value}.",
                 handle_cmdparam, (intptr_t)&cmdps, 0),
      OPT_INTEGER('t', "threads", &nr_threads,
                  "The number of threads to use on each MPI rank. Defaults to "
                  "1 if not specified.",
                  NULL, 0, 0),
      OPT_INTEGER('v', "verbose", &verbose,
                  "Run in verbose mode, in MPI mode 2 outputs from all ranks.",
                  NULL, 0, 0),
      OPT_INTEGER('y', "task-dumps", &dump_tasks,
                  "Time-step frequency at which task graphs are dumped.", NULL,
                  0, 0),
      OPT_INTEGER('Y', "threadpool-dumps", &dump_threadpool,
                  "Time-step frequency at which threadpool tasks are dumped.",
                  NULL, 0, 0),
      OPT_END(),
  };
  struct argparse argparse;
  argparse_init(&argparse, options, fof_usage, 0);
  argparse_describe(&argparse, "\nParameters:",
                    "\nSee the file examples/parameter_example.yml for an "
                    "example of parameter file.");
  int nargs = argparse_parse(&argparse, argc, (const char **)argv);

  /* Need a parameter file. */
  if (nargs != 1) {
    if (myrank == 0) argparse_usage(&argparse);
    pretime_message("\nError: no parameter file was supplied.\n");
    return 1;
  }
  param_filename = argv[0];

  /* Checks of options. */
#if !defined(HAVE_SETAFFINITY) || !defined(HAVE_LIBNUMA)
  if (with_aff) {
    pretime_message("Error: no NUMA support for thread affinity\n");
    return 1;
  }
#endif

#ifndef HAVE_FE_ENABLE_EXCEPT
  if (with_fp_exceptions) {
    pretime_message("Error: no support for floating point exceptions\n");
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

#ifndef SWIFT_DEBUG_THREADPOOL
  if (dump_threadpool) {
    pretime_message(
        "Error: threadpool dumping is only possible if SWIFT was "
        "configured with the --enable-threadpool-debugging option.\n");
    return 1;
  }
#endif

  /* The CPU frequency is a long long, so we need to parse that ourselves. */
  if (cpufreqarg != NULL) {
    if (sscanf(cpufreqarg, "%llu", &cpufreq) != 1) {
      if (myrank == 0)
        pretime_message("Error parsing CPU frequency (%s).\n", cpufreqarg);
      return 1;
    }
  }

  /* Write output parameter file */
  if (myrank == 0 && output_parameters_filename != NULL) {
    io_write_output_field_parameter(output_parameters_filename, with_cosmology,
                                    /*with_fof=*/1, /*with_stf=*/0);
    pretime_message("End of run.\n");
    return 0;
  }

/* Let's pin the main thread, now we know if affinity will be used. */
#if defined(HAVE_SETAFFINITY) && defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
  if (with_aff &&
      ((ENGINE_POLICY)&engine_policy_setaffinity) == engine_policy_setaffinity)
    engine_pin();
#endif

  /* Genesis 1.1: And then, there was time ! */
  clocks_set_cpufreq(cpufreq);

  /* How vocal are we ? */
  const int talking = (verbose == 1 && myrank == 0) || (verbose == 2);

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
    message("sizeof(gpart)       is %4zi bytes.", sizeof(struct gpart));
    message("sizeof(multipole)   is %4zi bytes.", sizeof(struct multipole));
    message("sizeof(grav_tensor) is %4zi bytes.", sizeof(struct grav_tensor));
    message("sizeof(task)        is %4zi bytes.", sizeof(struct task));
    message("sizeof(cell)        is %4zi bytes.", sizeof(struct cell));
  }

  /* Read the parameter file. */
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

  /* "Read" the Select Output file - this should actually do nothing, but
   * we need to mock the struct up for passing to `engine_init` */

  struct output_options *output_options =
      (struct output_options *)malloc(sizeof(struct output_options));
  output_options_init(params, myrank, output_options);

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

  /* Prepare and verify the selection of outputs */
  io_prepare_output_fields(output_options, with_cosmology, /*with_fof=*/1,
                           /*with_structure_finding=*/0, talking);

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
  const int clean_smoothing_length_values = parser_get_opt_param_int(
      params, "InitialConditions:cleanup_smoothing_lengths", 0);
  const int cleanup_h = parser_get_opt_param_int(
      params, "InitialConditions:cleanup_h_factors", 0);
  const int cleanup_sqrt_a = parser_get_opt_param_int(
      params, "InitialConditions:cleanup_velocity_factors", 0);

  /* Initialise the cosmology */
  if (with_cosmology)
    cosmology_init(params, &us, &prog_const, &cosmo);
  else
    cosmology_init_no_cosmo(&cosmo);
  if (myrank == 0 && with_cosmology) cosmology_print(&cosmo);

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

  /* Initialise the FOF properties */
  bzero(&fof_properties, sizeof(struct fof_props));
  fof_init(&fof_properties, params, &prog_const, &us, /*stand-alone=*/1);

  /* Be verbose about what happens next */
  if (myrank == 0) message("Reading ICs from file '%s'", ICfileName);
  if (myrank == 0 && cleanup_h)
    message("Cleaning up h-factors (h=%f)", cosmo.h);
  if (myrank == 0 && cleanup_sqrt_a)
    message("Cleaning up a-factors from velocity (a=%f)", cosmo.a);
  fflush(stdout);

  /* Get ready to read particles of all kinds */
  int flag_entropy_ICs = 0;
  size_t Ngas = 0, Ngpart = 0, Ngpart_background = 0, Nnupart = 0;
  size_t Nsink = 0, Nspart = 0, Nbpart = 0;
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
                   /*with_grav=*/1, with_sinks, with_stars, with_black_holes,
                   with_cosmology, cleanup_h, cleanup_sqrt_a, cosmo.h, cosmo.a,
                   myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL, nr_threads,
                   /*dry_run=*/0, /*remap_ids=*/0, &ics_metadata);
#else
  read_ic_serial(ICfileName, &us, dim, &parts, &gparts, &sinks, &sparts,
                 &bparts, &Ngas, &Ngpart, &Ngpart_background, &Nnupart, &Nsink,
                 &Nspart, &Nbpart, &flag_entropy_ICs, with_hydro,
                 /*with_grav=*/1, with_sinks, with_stars, with_black_holes,
                 with_cosmology, cleanup_h, cleanup_sqrt_a, cosmo.h, cosmo.a,
                 myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL, nr_threads,
                 /*dry_run=*/0, /*remap_ids=*/0, &ics_metadata);
#endif
#else
  read_ic_single(ICfileName, &us, dim, &parts, &gparts, &sinks, &sparts,
                 &bparts, &Ngas, &Ngpart, &Ngpart_background, &Nnupart, &Nsink,
                 &Nspart, &Nbpart, &flag_entropy_ICs, with_hydro,
                 /*with_grav=*/1, with_sinks, with_stars, with_black_holes,
                 with_cosmology, cleanup_h, cleanup_sqrt_a, cosmo.h, cosmo.a,
                 nr_threads, /*dry_run=*/0, /*remap_ids=*/0, &ics_metadata);
#endif
#endif
  if (myrank == 0) {
    clocks_gettime(&toc);
    message("Reading initial conditions took %.3f %s.", clocks_diff(&tic, &toc),
            clocks_getunit());
    fflush(stdout);
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that the other links are correctly set */
  part_verify_links(parts, gparts, sinks, sparts, bparts, Ngas, Ngpart, Nsink,
                    Nspart, Nbpart, /*verbose=*/1);
#endif

  /* Get the total number of particles across all nodes. */
  long long N_total[swift_type_count + 1] = {0};
  long long Nbaryons = Ngas + Nspart + Nbpart;
#if defined(WITH_MPI)
  long long N_long[swift_type_count + 1] = {0};
  N_long[swift_type_gas] = Ngas;
  N_long[swift_type_dark_matter] =
      Ngpart - Ngpart_background - Nnupart - Nbaryons;
  N_long[swift_type_dark_matter_background] = Ngpart_background;
  N_long[swift_type_sink] = Nsink;
  N_long[swift_type_stars] = Nspart;
  N_long[swift_type_black_hole] = Nbpart;
  N_long[swift_type_neutrino] = Nnupart;
  N_long[swift_type_count] = Ngpart;
  MPI_Allreduce(N_long, N_total, swift_type_count + 1, MPI_LONG_LONG_INT,
                MPI_SUM, MPI_COMM_WORLD);
#else
  N_total[swift_type_gas] = Ngas;
  N_total[swift_type_dark_matter] =
      Ngpart - Ngpart_background - Nnupart - Nbaryons;
  N_total[swift_type_dark_matter_background] = Ngpart_background;
  N_total[swift_type_sink] = Nsink;
  N_total[swift_type_stars] = Nspart;
  N_total[swift_type_black_hole] = Nbpart;
  N_total[swift_type_neutrino] = Nnupart;
  N_total[swift_type_count] = Ngpart;
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
      (N_total[swift_type_gas] + N_total[swift_type_sink] +
       N_total[swift_type_stars] + N_total[swift_type_black_hole]) > 0;

  /* Do we have background DM particles? */
  const int with_DM_background_particles =
      N_total[swift_type_dark_matter_background] > 0;

  /* Do we have neutrino DM particles? */
  const int with_neutrinos = N_total[swift_type_neutrino] > 0;

  /* Zero out neutrino properties to avoid running neutrino tasks */
  bzero(&neutrino_properties, sizeof(struct neutrino_props));

  /* Initialize the space with these data. */
  if (myrank == 0) clocks_gettime(&tic);
  space_init(&s, params, &cosmo, dim, /*hydro_props=*/NULL, parts, gparts,
             sinks, sparts, bparts, Ngas, Ngpart, Nsink, Nspart, Nbpart,
             Nnupart, periodic, replicate, /*remap_ids=*/0,
             /*generate_gas_in_ics=*/0, /*hydro=*/N_total[0] > 0, /*gravity=*/1,
             /*with_star_formation=*/0, /*sink=*/N_total[swift_type_sink],
             with_DM_particles, with_DM_background_particles, with_neutrinos,
             talking, /*dry_run=*/0, nr_nodes);

  if (myrank == 0) {
    clocks_gettime(&toc);
    message("space_init took %.3f %s.", clocks_diff(&tic, &toc),
            clocks_getunit());
    fflush(stdout);
  }

  /* Initialise the gravity scheme */
  bzero(&gravity_properties, sizeof(struct gravity_props));
  gravity_props_init(&gravity_properties, params, &prog_const, &cosmo,
                     with_cosmology, /*with_external_gravity=*/0,
                     with_baryon_particles, with_DM_particles, with_neutrinos,
                     with_DM_background_particles, periodic, s.dim, s.cdim);

  /* Initialise the long-range gravity mesh */
  if (periodic) {
#ifdef HAVE_FFTW
    pm_mesh_init(&mesh, &gravity_properties, s.dim, nr_threads);
#else
    /* Need the FFTW library if periodic and self gravity. */
    error("No FFTW library found. Cannot compute periodic long-range forces.");
#endif
  } else {
    pm_mesh_init_no_mesh(&mesh, s.dim);
  }

  /* Also update the total counts (in case of changes due to replication) */
  Nbaryons = s.nr_parts + s.nr_sparts + s.nr_bparts;
  Nnupart = s.nr_nuparts;
#if defined(WITH_MPI)
  N_long[swift_type_gas] = s.nr_parts;
  N_long[swift_type_dark_matter] =
      s.nr_gparts - Ngpart_background - Nbaryons - Nnupart;
  N_long[swift_type_count] = s.nr_gparts;
  N_long[swift_type_sink] = s.nr_sinks;
  N_long[swift_type_stars] = s.nr_sparts;
  N_long[swift_type_black_hole] = s.nr_bparts;
  N_long[swift_type_neutrino] = s.nr_nuparts;
  MPI_Allreduce(N_long, N_total, swift_type_count + 1, MPI_LONG_LONG_INT,
                MPI_SUM, MPI_COMM_WORLD);
#else
  N_total[swift_type_gas] = s.nr_parts;
  N_total[swift_type_dark_matter] =
      s.nr_gparts - Ngpart_background - Nbaryons - Nnupart;
  N_total[swift_type_count] = s.nr_gparts;
  N_total[swift_type_sink] = s.nr_sinks;
  N_total[swift_type_stars] = s.nr_sparts;
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
    message("maximum depth is %d.", s.maxdepth);
    fflush(stdout);
  }

  /* Construct the engine policy */
  int engine_policies = ENGINE_POLICY | engine_policy_steal;
  engine_policies |= engine_policy_self_gravity;
  engine_policies |= engine_policy_fof;
  if (with_cosmology) engine_policies |= engine_policy_cosmology;

  /* Initialize the engine with the space and policies. */
  engine_init(
      &e, &s, params, output_options, N_total[swift_type_gas],
      N_total[swift_type_count], N_total[swift_type_sink],
      N_total[swift_type_stars], N_total[swift_type_black_hole],
      N_total[swift_type_dark_matter_background], N_total[swift_type_neutrino],
      engine_policies, talking, &us, &prog_const, &cosmo, &hydro_properties,
      /*entropy_floor=*/NULL, &gravity_properties,
      /*stars_properties=*/NULL, /*black_holes_properties=*/NULL,
      /*sink_properties=*/NULL, &neutrino_properties,
      /*neutrino_response=*/NULL, /*feedback_properties=*/NULL,
      /*pressure_floor_properties=*/NULL,
      /*rt_properties=*/NULL, &mesh, /*pow_data=*/NULL, /*potential=*/NULL,
      /*forcing_terms=*/NULL,
      /*cooling_func=*/NULL, /*starform=*/NULL, /*chemistry=*/NULL,
      /*extra_io_props=*/NULL, &fof_properties, /*los_properties=*/NULL,
      /*lightcone_properties=*/NULL, &ics_metadata);
  engine_config(/*restart=*/0, /*fof=*/1, &e, params, nr_nodes, myrank,
                nr_threads, nr_threads, with_aff, talking, NULL, NULL,
                &reparttype);

  /* Get some info to the user. */
  if (myrank == 0) {
    const long long N_DM = N_total[swift_type_dark_matter] +
                           N_total[swift_type_dark_matter_background];

    message(
        "Running FOF on %lld gas particles, %lld sink particles, %lld stars "
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

#ifdef WITH_MPI
  /* Split the space. */
  engine_split(&e, &initial_partition);
#endif

#ifdef SWIFT_DEBUG_TASKS
  e.tic_step = getticks();
#endif

  /* Initialise the tree and communication tasks */
  engine_rebuild(&e, /*repartitionned=*/0, clean_smoothing_length_values);

#ifdef SWIFT_DEBUG_TASKS
  e.toc_step = getticks();
#endif

  /* Perform the FOF search */
  engine_fof(&e, /*dump_results=*/1, /*dump_debug=*/0, /*seed_black_holes=*/0,
             /*buffers allocated=*/1);

  /* Update the policies to make sure the particles are written
   * if they exist */
  if (with_hydro) e.policy |= engine_policy_hydro;
  if (with_stars) e.policy |= engine_policy_stars;
  if (with_black_holes) e.policy |= engine_policy_black_holes;
  if (with_sinks) e.policy |= engine_policy_sinks;

  /* Write output. */
  engine_dump_snapshot(&e, /*fof=*/1);

#ifdef WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* Dump the task data using the given frequency. */
  if (dump_tasks) {
#ifdef SWIFT_DEBUG_TASKS
    task_dump_all(&e, 0);
#endif

    /* Generate the task statistics. */
    char dumpfile[40];
    snprintf(dumpfile, 40, "thread_stats-step%d.dat", 0);
    task_dump_stats(dumpfile, &e, /* dump_tasks_threshold = */ 0.f,
                    /* header = */ 0, /* allranks = */ 1);
  }

#ifdef SWIFT_DEBUG_THREADPOOL
  /* Dump the task data using the given frequency. */
  if (dump_threadpool) {
    char threadpool_dumpfile[52];
#ifdef WITH_MPI
    snprintf(threadpool_dumpfile, 52, "threadpool_info-rank%d-step%d.dat",
             engine_rank, 0);
#else
    snprintf(threadpool_dumpfile, 52, "threadpool_info-step%d.dat", 0);
#endif  // WITH_MPI
    threadpool_dump_log(&e.threadpool, threadpool_dumpfile, 1);
  } else {
    threadpool_reset_log(&e.threadpool);
  }
#endif  // SWIFT_DEBUG_THREADPOOL

  /* used parameters */
  parser_write_params_to_file(params, "fof_used_parameters.yml", /*used=*/1);
  /* unused parameters */
  parser_write_params_to_file(params, "fof_unused_parameters.yml", /*used=*/0);

  /* Dump memory use report */
#ifdef SWIFT_MEMUSE_REPORTS
  {
    char dumpfile[40];
#ifdef WITH_MPI
    snprintf(dumpfile, 40, "memuse_report-rank%d-fof.dat", engine_rank);
#else
    snprintf(dumpfile, 40, "memuse_report-fof.dat");
#endif  // WITH_MPI
    memuse_log_dump(dumpfile);
  }
#endif

  /* Clean everything */
  cosmology_clean(&cosmo);
  pm_mesh_clean(&mesh);
  engine_clean(&e, /*fof=*/1, /*restart=*/0);
  free(params);
  free(output_options);

#ifdef WITH_MPI
  partition_clean(&initial_partition, &reparttype);
  if ((res = MPI_Finalize()) != MPI_SUCCESS)
    error("call to MPI_Finalize failed with error %i.", res);
#endif

  /* Say goodbye. */
  if (myrank == 0) message("done. Bye.");

  /* All is calm, all is bright. */
  return 0;
}
