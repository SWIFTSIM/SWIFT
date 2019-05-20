/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 STFC (author: aidan.chalk@stfc.ac.uk)
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

//TODO Option for parameter overrides
/*static int handle_cmdparam(struct argparse *self,
                           const struct argparse_option *opt) {
  struct cmdparams *cmdps = (struct cmdparams *)opt->data;
  cmdps->param[cmdps->nparam] = *(char **)opt->value;
  cmdps->nparam++;
  return 1;
}*/



    

int main(int argc, char *argv[]){

struct gpart *gparts = NULL;
struct part *parts = NULL;
struct spart *sparts = NULL;
struct bpart *bparts = NULL;
struct unit_system us;
struct space s;
struct engine e;

int flag_entropy_ICs = 0;
const int cleanup_h = 0;
const int cleanup_sqrt_a = 0;

size_t Ngas = 0, Ngpart = 0, Nspart = 0, Nboundary = 0, Nfluid = 0, Nblackhole=0;
double dim[3] = {0., 0., 0.};
int nr_nodes = 1, myrank = 0;

#if defined(WITH_MPI)
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

units_init_cgs(&us);

/* Read command line options */
int with_aff = 0;
int dump_tasks = 0;
int dump_threadpool = 0;
int with_fp_exceptions = 0;
char *cpufreqarg = NULL;
int nr_threads = 1;
int verbose = 0;
struct cmdparams cmdps;
cmdps.nparam = 0;
cmdps.param[0] = NULL;
unsigned long long cpufreq = 0;
char *param_filename = NULL;
struct phys_const prog_const;
struct pm_mesh mesh;
char restart_file[200] = "";
int clean_smoothing_length_values = 0;



struct argparse_option options[] = {
    OPT_HELP(),
    OPT_GROUP("   Control options:\n"),
    OPT_BOOLEAN('a', "pin", &with_aff,
                "Pin runners using processor affinity.", NULL, 0, 0),
    OPT_BOOLEAN('e', "fpe", &with_fp_exceptions,
                "Enable floating-point exceptions (debugging mode).", NULL, 0,
                 0),
    OPT_STRING('f', "cpu-frequency", &cpufreqarg,
               "Overwrite the CPU "
               "frequency (Hz) to be used for time measurements.",
                NULL, 0, 0),
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
argparse_init(&argparse, options, swift_usage, 0);
argparse_describe(&argparse, "\nParameters:",
                  "\nSee the file examples/parameter_example.yml for an "
                  "example of parameter file.");
int nargs = argparse_parse(&argparse, argc, (const char **)argv);

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

#ifndef HAVE_FE_ENABLE_EXCEPT
  if (with_fp_exceptions) {
    printf("Error: no support for floating point exceptions\n");
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
    printf(
        "Error: threadpool dumping is only possible if SWIFT was "
        "configured with the --enable-threadpool-debugging option.\n");
    return 1;
  }
#endif

  /* The CPU frequency is a long long, so we need to parse that ourselves. */
  if (cpufreqarg != NULL) {
    if (sscanf(cpufreqarg, "%llu", &cpufreq) != 1) {
      if (myrank == 0)
        printf("Error parsing CPU frequency (%s).\n", cpufreqarg);
      return 1;
    }
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
    message("sizeof(spart)       is %4zi bytes.", sizeof(struct spart));
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

  /* Check that we can write the snapshots by testing if the output
 *    * directory exists and is searchable and writable. */
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
    const int generate_gas_in_ics = 0;

   /* No cosmology with engineering */
   struct cosmology cosmo; 
   cosmology_init_no_cosmo(&cosmo);

   /* No stars with engineering */
   struct stars_props stars_properties;
   bzero(&stars_properties, sizeof(struct stars_props));

   /* No gravity with engineering */
   struct gravity_props gravity_properties;
   bzero(&gravity_properties, sizeof(struct gravity_props));

   /* No external potential properties with engineering */
   struct external_potential potential;
   bzero(&potential, sizeof(struct external_potential));

   /*No cosmological cooling with engineering */
   struct cooling_function_data cooling_func;
   bzero(&cooling_func, sizeof(struct cooling_function_data));
 
   /* No star formation with engineering */
   struct star_formation starform;
   bzero(&starform, sizeof(struct star_formation));

   /* No cosmological chemistry with engineering */
   struct chemistry_global_data chemistry;
   bzero(&chemistry, sizeof(struct chemistry_global_data));

   /* Nullify hydro properties for now - these need changing for Engineering SPH */
   struct hydro_props hydro_properties;
   hydro_props_init(&hydro_properties, &prog_const, &us, params);

   /* Nullify equation of state for now - these need changing for Engineering SPH */
   eos_init(&eos, &prog_const, &us, params);

   /* Nullify entropy floor for now, no idea what it is. */
   struct entropy_floor_properties entropy_floor;
   bzero(&entropy_floor, sizeof(struct entropy_floor_properties));

  struct black_holes_props black_hole_properties;
  bzero(&black_hole_properties, sizeof(struct black_holes_props));


  struct feedback_props feedback_properties;
  bzero(&feedback_properties, sizeof(struct feedback_props));

  /* Be verbose about what happens next */
  if (myrank == 0) message("Reading ICs from file '%s'", ICfileName);
  fflush(stdout);


#if defined(HAVE_HDF5)
#if defined(WITH_MPI)
    read_ic_serial("boundary_test.hdf5",&us, dim, &parts, &gparts, &sparts, &bparts, &Ngas,
                   &Ngpart, &Nspart, &Nboundary, &Nblackhole, &Nfluid, &flag_entropy_ICs, 0,
                   0, 0, 0, 1, cleanup_h, cleanup_sqrt_a, 1.0, 1.0, myrank, nr_nodes, MPI_COMM_WORLD,
                   MPI_INFO_NULL, 1, 0);
#else
    read_ic_single(ICfileName,&us, dim, &parts, &gparts, &sparts, &bparts, &Ngas,
                   &Ngpart, &Nspart, &Nblackhole, &Nboundary, &Nfluid, &flag_entropy_ICs, 0,
                   0, 0, 0, 1, cleanup_h, cleanup_sqrt_a, 1.0, 1.0, 1, 0);
#endif
#else
    error("Failed to find MPI and HDF5");
#endif

long long N_total[3] = {0, 0, 0};
#if defined(WITH_MPI)
    long long N_long[3] = {Ngas, Ngpart, Nspart};
    MPI_Allreduce(&N_long, &N_total, 3, MPI_LONG_LONG_INT, MPI_SUM,
                  MPI_COMM_WORLD);
#else
    N_total[0] = Nfluid + Nboundary;
    N_total[1] = 0;
    N_total[2] = 0;
#endif

/* Engineering doesn't read gas, dark matter or star particles, so lets check. */
if(Ngas != 0){
    error("Found gas particles");
}
if(Ngpart != 0){
    error("Found gravity particles");
}
if(Nspart != 0){
    error("Found stars");
}

if(Nblackhole != 0){
    error("Found blackholes");
}

if (myrank == 0){
  message( "Read %lu boundary particles and a total of %lli particles", Nboundary, N_total[0] );
}

//TODO Remove these checks.
/*if(Nboundary != 6){
    error("Found %lu boundary particles", Nboundary);
}
if(Nfluid != 1){
    error("Found %lu fluid particles", Nfluid);
}

if(Nboundary == 6 && Nfluid  == 1){
    message("We found the right number of particles!");
}*/

/*for(size_t i = 0; i < Nboundary; i++){
    printf("%f %f %f %f %i %i\n", parts[i].x[1], parts[i].v[1], parts[i].h, parts[i].rho, part_is_active(&parts[i], &e), parts[i].is_boundary);
}
for(size_t i = Nboundary; i < Nboundary+Nfluid; i++){
    printf("%f %f %f %f %i %i\n", parts[i].x[1], parts[i].v[1], parts[i].h, parts[i].rho, part_is_active(&parts[i], &e), parts[i].is_boundary);
}*/

//TODO Parse parameter
/* (myrank == 0){
  io_check_output_fields(params, N_total);
}*/

space_init(&s, params, &cosmo, dim, parts, gparts, sparts, bparts, Nfluid+Nboundary, Ngpart,
           Nspart, Nblackhole, periodic, replicate, generate_gas_in_ics, 1 /*with_hydro*/, 0 /*with_self_gravity*/,
           0 /*with_star_formation*/, talking, 0 /*dry_run*/); 


   /* No gravity with engineering */
   pm_mesh_init_no_mesh(&mesh, s.dim);



int engine_policies = ENGINE_POLICY | engine_policy_steal;
engine_policies |= engine_policy_hydro;
engine_policies |= engine_policy_engineering;



engine_init(&e, &s, params, N_total[0], N_total[1], N_total[2],
            engine_policies, talking, &reparttype, &us, &prog_const, &cosmo,
            &hydro_properties, &entropy_floor, &gravity_properties,
            &stars_properties, &black_hole_properties, &feedback_properties, &mesh, &potential, &cooling_func, &starform,
            &chemistry);
engine_config(0, &e, params, nr_nodes, myrank, nr_threads, with_aff,
              talking, restart_file);

#ifdef WITH_MPI
/* Split the space. */
engine_split(&e, &initial_partition);
engine_redistribute(&e);
#endif

/* Initialise the particles */
engine_init_particles(&e, flag_entropy_ICs, clean_smoothing_length_values);

/* Write the state of the system before starting time integration. */
engine_dump_snapshot(&e);
engine_print_stats(&e);

  /* Legend */
  if (myrank == 0) {
    printf("# %6s %14s %12s %12s %14s %9s %12s %12s %12s %12s %16s [%s] %6s\n",
           "Step", "Time", "Scale-factor", "Redshift", "Time-step", "Time-bins",
           "Updates", "g-Updates", "s-Updates", "b-Updates", "Wall-clock time",
           clocks_getunit(), "Props");
    fflush(stdout);
  }

  /* dump the parameters as used. */

  /* used parameters */
  parser_write_params_to_file(params, "used_parameters.yml", 1);
  /* unused parameters */
  parser_write_params_to_file(params, "unused_parameters.yml", 0);
#ifdef HAVE_FE_ENABLE_EXCEPT
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  for (int j = 0; !engine_is_done(&e) && e.step -1 != 2202199; j++){
    timers_reset_all();

    engine_step(&e);
for(size_t i = 0; i < Nboundary+Nfluid; i++){
  if(!parts[i].is_boundary){
   printf("%f %f\n", parts[i].a_hydro[0], parts[i].a_constant[0]);
    break;
  }
}
for(size_t i = 0; i < Nboundary; i++){
    //printf("%f %f %f %f %f %i\n",parts[i].x[0], parts[i].x[1], /*parts[i].drho_dt,*/ parts[i].v[1], parts[i].h, parts[i].rho, part_is_active(&parts[i], &e)/*, parts[i].is_boundary*/);
  //  printf("%f %f %f %f %f %f %i %i\n",parts[i].x[0], parts[i].x[1], parts[i].drho_dt, parts[i].v[1], parts[i].h, parts[i].rho, part_is_active(&parts[i], &e), parts[i].is_boundary);
}
for(size_t i = Nboundary; i < Nboundary+Nfluid; i++){
    //printf("%f %f %f %f %f %i\n",parts[i].x[0], parts[i].x[1], /*parts[i].drho_dt,*/ parts[i].v[1], parts[i].h, parts[i].rho, part_is_active(&parts[i], &e)/*, parts[i].is_boundary*/);
//    printf("%f %f %f %f %f %f %i %i\n",parts[i].x[0], parts[i].x[1], parts[i].drho_dt, parts[i].v[1], parts[i].h, parts[i].rho, part_is_active(&parts[i], &e), parts[i].is_boundary);
}

/*for(int i = 0; i < e.s->nr_cells_with_particles; i++){
  struct cell *c = &e.s->cells_top[e.s->cells_with_particles_top[i]];
  printf("cell %f %f %f has %i particles\n", c->loc[0], c->loc[1], c->loc[2], c->hydro.count);
  if(c->loc[0] == 30.0 && c->loc[1] == 90.0 && c->loc[2] == 30.0){
    struct link *z = c->hydro.density;
    int count = 0;
    while( z != NULL ){
      count++;
      z = z->next;
    }
  printf("cell we're looking for has %i density tasks\n", count);
  printf("Does cell violate rule: %f %f\n", kernel_gamma * c->hydro.h_max + c->hydro.dx_max_part, c->dmin);
  }
}*/

//    engine_dump_snapshot(&e);
  }

    engine_drift_all(&e, /*drift_mpole=*/0);
    engine_print_stats(&e);
    /* write a final snapshot */
    engine_dump_snapshot(&e);

printf("engine is done %i e.step %i\n", engine_is_done(&e), e.step);
printf("e.ti_current = %lli\n", e.ti_current);

/*TODO Work around compile error*/
message("%i %i %i", nr_nodes, talking, clean_smoothing_length_values);

#if defined(WITH_MPI)
MPI_Finalize();
#endif
  engine_clean(&e);
free(params);


 /* Say goodbye. */
  if (myrank == 0) message("done. Bye.");


return 0;



}
