/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "../config.h"

/* Some standard headers. */
#include <fenv.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

/**
 * @brief Help messages for the command line parameters.
 */
void print_help_message() {

  printf("\nUsage: swift [OPTION] PARAMFILE\n\n");
  printf("Valid options are:\n");
  printf("  %2s %8s %s\n", "-c", "", "Run with cosmological time integration");
  printf("  %2s %8s %s\n", "-e", "",
         "Enable floating-point exceptions (debugging mode)");
  printf("  %2s %8s %s\n", "-f", "[value] ",
         "Overwrites the CPU frequency (Hz) to be used for time measurements");
  printf("  %2s %8s %s\n", "-g", "",
         "Run with an external gravitational potential");
  printf("  %2s %8s %s\n", "-G", "", "Run with self-gravity");
  printf("  %2s %8s %s\n", "-s", "", "Run with SPH");
  printf("  %2s %8s %s\n", "-v", "[12]   ",
         "Increase the level of verbosity 1: MPI-rank 0 writes "
         "2: All MPI ranks write");
  printf("  %2s %8s %s\n", "-y", "",
         "Time-step frequency at which task graphs are dumped");
  printf("  %2s %8s %s\n", "-h", "", "Prints this help message and exits");
  printf("\nSee the file example.yml for an example of parameter file.\n");
}

/**
 * @brief Main routine that loads a few particles and generates some output.
 *
 */
int main(int argc, char *argv[]) {

  struct clocks_time tic, toc;

#ifdef WITH_MPI
  struct partition initial_partition;
  enum repartition_type reparttype = REPART_NONE;

  initial_partition.type = INITPART_GRID;
  initial_partition.grid[0] = 1;
  initial_partition.grid[1] = 1;
  initial_partition.grid[2] = 1;
#ifdef HAVE_METIS
  /* Defaults make use of METIS. */
  reparttype = REPART_METIS_BOTH;
  initial_partition.type = INITPART_METIS_NOWEIGHT;
#endif
#endif

  int nr_nodes = 1, myrank = 0;
#ifdef WITH_MPI
  /* Start by initializing MPI. */
  int res = 0, prov = 0;
  if ((res = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &prov)) !=
      MPI_SUCCESS)
    error("Call to MPI_Init failed with error %i.", res);
  if (prov != MPI_THREAD_MULTIPLE)
    error(
        "MPI does not provide the level of threading required "
        "(MPI_THREAD_MULTIPLE).");
  if ((res = MPI_Comm_size(MPI_COMM_WORLD, &nr_nodes)) != MPI_SUCCESS)
    error("MPI_Comm_size failed with error %i.", res);
  if ((res = MPI_Comm_rank(MPI_COMM_WORLD, &myrank)) != MPI_SUCCESS)
    error("Call to MPI_Comm_rank failed with error %i.", res);
  if ((res = MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN)) !=
      MPI_SUCCESS)
    error("Call to MPI_Comm_set_errhandler failed with error %i.", res);
  if (myrank == 0)
    printf("[0000][00000.0] MPI is up and running with %i node(s).\n",
           nr_nodes);
  fflush(stdout);

  /* Set a default grid so that grid[0]*grid[1]*grid[2] == nr_nodes. */
  factor(nr_nodes, &initial_partition.grid[0], &initial_partition.grid[1]);
  factor(nr_nodes / initial_partition.grid[1], &initial_partition.grid[0],
         &initial_partition.grid[2]);
  factor(initial_partition.grid[0] * initial_partition.grid[1],
         &initial_partition.grid[1], &initial_partition.grid[0]);
#endif

  /* Greeting message */
  if (myrank == 0) greetings();

#if defined(HAVE_SETAFFINITY) && defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)
  if ((ENGINE_POLICY) & engine_policy_setaffinity) {
    /* Ensure the NUMA node on which we initialise (first touch) everything
     * doesn't change before engine_init allocates NUMA-local workers.
     * Otherwise, we may be scheduled elsewhere between the two times.
     */
    cpu_set_t affinity;
    CPU_ZERO(&affinity);
    CPU_SET(sched_getcpu(), &affinity);
    if (sched_setaffinity(0, sizeof(cpu_set_t), &affinity) != 0) {
      error("failed to set entry thread's affinity");
    }
  }
#endif

  int dump_tasks = 0;
  int with_cosmology = 0;
  int with_external_gravity = 0;
  int with_self_gravity = 0;
  int with_hydro = 0;
  int with_fp_exceptions = 0;
  int verbose = 0;
  char paramFileName[200] = "";
  unsigned long long cpufreq = 0;

  /* Parse the parameters */
  int c;
  while ((c = getopt(argc, argv, "cef:gGhsv:y")) != -1) switch (c) {
      case 'c':
        with_cosmology = 1;
        break;
      case 'e':
        with_fp_exceptions = 1;
        break;
      case 'f':
        if (sscanf(optarg, "%llu", &cpufreq) != 1) {
          if (myrank == 0) printf("Error parsing CPU frequency (-f).\n");
          if (myrank == 0) print_help_message();
          exit(1);
        }
        break;
      case 'g':
        with_external_gravity = 1;
        break;
      case 'G':
        with_self_gravity = 1;
        break;
      case 'h':
        if (myrank == 0) print_help_message();
        exit(0);
      case 's':
        with_hydro = 1;
        break;
      case 'v':
        if (sscanf(optarg, "%d", &verbose) != 1) {
          if (myrank == 0) printf("Error parsing verbosity level (-v).\n");
          if (myrank == 0) print_help_message();
          exit(1);
        }
        break;
      case 'y':
        if (sscanf(optarg, "%d", &dump_tasks) != 1) {
          if (myrank == 0) printf("Error parsing dump_tasks (-y). \n");
          if (myrank == 0) print_help_message();
          exit(1);
        }
        break;
      case '?':
        if (myrank == 0) print_help_message();
        exit(1);
        break;
    }
  if (optind == argc - 1) {
    if (!strcpy(paramFileName, argv[optind++]))
      error("Error reading parameter file name.");
  } else if (optind > argc - 1) {
    if (myrank == 0) printf("Error: A parameter file name must be provided\n");
    if (myrank == 0) print_help_message();
    exit(1);
  } else {
    if (myrank == 0) printf("Error: Too many parameters given\n");
    if (myrank == 0) print_help_message();
    exit(1);
  }

  /* And then, there was time ! */
  clocks_set_cpufreq(cpufreq);

  /* Report CPU frequency. */
  if (myrank == 0) {
    cpufreq = clocks_get_cpufreq();
    message("CPU frequency used for tick conversion: %llu Hz", cpufreq);
  }

  /* Do we choke on FP-exceptions ? */
  if (with_fp_exceptions) {
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    if (myrank == 0) message("Floating point exceptions will be reported.");
  }

  /* How large are the parts? */
  if (myrank == 0) {
    message("sizeof(struct part)  is %4zi bytes.", sizeof(struct part));
    message("sizeof(struct xpart) is %4zi bytes.", sizeof(struct xpart));
    message("sizeof(struct gpart) is %4zi bytes.", sizeof(struct gpart));
  }

  /* How vocal are we ? */
  const int talking = (verbose == 1 && myrank == 0) || (verbose == 2);

  /* Read the parameter file */
  struct swift_params params;
  if (myrank == 0) {
    message("Reading parameters from file '%s'", paramFileName);
    parser_read_file(paramFileName, &params);
    // parser_print_params(&params);
    parser_write_params_to_file(&params, "used_parameters.yml");
  }
#ifdef WITH_MPI
  /* Broadcast the parameter file */
  MPI_Bcast(&params, sizeof(struct swift_params), MPI_BYTE, MPI_COMM_WORLD);
#endif

  /* Initialize unit system */
  struct UnitSystem us;
  units_init(&us, &params);
  if (myrank == 0) {
    message("Unit system: U_M = %e g.", us.UnitMass_in_cgs);
    message("Unit system: U_L = %e cm.", us.UnitLength_in_cgs);
    message("Unit system: U_t = %e s.", us.UnitTime_in_cgs);
    message("Unit system: U_I = %e A.", us.UnitCurrent_in_cgs);
    message("Unit system: U_T = %e K.", us.UnitTemperature_in_cgs);
  }

/* Some initial information about domain decomposition */
#ifdef WITH_MPI
  if (myrank == 0) {
    // message("Running with %i thread(s) per node.", nr_threads);
    message("Using initial partition %s",
            initial_partition_name[initial_partition.type]);
    if (initial_partition.type == INITPART_GRID)
      message("grid set to [ %i %i %i ].", initial_partition.grid[0],
              initial_partition.grid[1], initial_partition.grid[2]);
    message("Using %s repartitioning", repartition_name[reparttype]);

    if (nr_nodes == 1) {
      message("WARNING: you are running with one MPI rank.");
      message("WARNING: you should use the non-MPI version of this program.");
    }
    fflush(stdout);
  }
#else
// if (myrank == 0) message("Running with %i thread(s).", nr_threads);
#endif

  /* Read particles and space information from (GADGET) ICs */
  char ICfileName[200] = "";
  parser_get_param_string(&params, "InitialConditions:file_name", ICfileName);
  struct part *parts = NULL;
  struct gpart *gparts = NULL;
  size_t Ngas = 0, Ngpart = 0;
  double dim[3] = {0., 0., 0.};
  int periodic = 0;
  if (myrank == 0) clocks_gettime(&tic);
#if defined(WITH_MPI)
#if defined(HAVE_PARALLEL_HDF5)
  read_ic_parallel(ICfileName, dim, &parts, &gparts, &Ngas, &Ngpart, &periodic,
                   myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
  read_ic_serial(ICfileName, dim, &parts, &gparts, &Ngas, &Ngpart, &periodic,
                 myrank, nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
#else
  read_ic_single(ICfileName, dim, &parts, &gparts, &Ngas, &Ngpart, &periodic);
#endif
  if (myrank == 0) {
    clocks_gettime(&toc);
    message("Reading initial conditions took %.3f %s.", clocks_diff(&tic, &toc),
            clocks_getunit());
    fflush(stdout);
  }

  /* Discard gparts if we don't have gravity
   * (Better implementation of i/o will come)*/
  if (!with_external_gravity && !with_self_gravity) {
    free(gparts);
    gparts = NULL;
    for (size_t k = 0; k < Ngas; ++k) parts[k].gpart = NULL;
    Ngpart = 0;
  }

  /* Get the total number of particles across all nodes. */
  long long N_total[2] = {0, 0};
#if defined(WITH_MPI)
  long long N_long[2] = {Ngas, Ngpart};
  MPI_Reduce(&N_long, &N_total, 2, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  N_total[0] = Ngas;
  N_total[1] = Ngpart;
#endif
  if (myrank == 0)
    message("Read %lld gas particles and %lld gparts from the ICs.", N_total[0],
            N_total[1]);

  /* Initialize the space with these data. */
  if (myrank == 0) clocks_gettime(&tic);
  struct space s;
  space_init(&s, &params, dim, parts, gparts, Ngas, Ngpart, periodic, talking);
  if (talking) {
    clocks_gettime(&toc);
    message("space_init took %.3f %s.", clocks_diff(&tic, &toc),
            clocks_getunit());
    fflush(stdout);
  }

  /* Say a few nice things about the space we just created. */
  if (talking) {
    message("space dimensions are [ %.3f %.3f %.3f ].", s.dim[0], s.dim[1],
            s.dim[2]);
    message("space %s periodic.", s.periodic ? "is" : "isn't");
    message("highest-level cell dimensions are [ %i %i %i ].", s.cdim[0],
            s.cdim[1], s.cdim[2]);
    message("%zi parts in %i cells.", s.nr_parts, s.tot_cells);
    message("%zi gparts in %i cells.", s.nr_gparts, s.tot_cells);
    message("maximum depth is %d.", s.maxdepth);
  }

  /* Verify that each particle is in it's proper cell. */
  if (talking) {
    int icount = 0;
    space_map_cells_pre(&s, 0, &map_cellcheck, &icount);
    message("map_cellcheck picked up %i parts.", icount);
  }

  /* Verify the maximal depth of cells. */
  if (talking) {
    int data[2] = {s.maxdepth, 0};
    space_map_cells_pre(&s, 0, &map_maxdepth, data);
    message("nr of cells at depth %i is %i.", data[0], data[1]);
  }

  /* Construct the engine policy */
  int engine_policies = ENGINE_POLICY | engine_policy_steal;
  if (with_hydro) engine_policies |= engine_policy_hydro;
  if (with_self_gravity) engine_policies |= engine_policy_self_gravity;
  if (with_external_gravity) engine_policies |= engine_policy_external_gravity;
  if (with_cosmology) engine_policies |= engine_policy_cosmology;

  /* Initialize the engine with the space and policies. */
  if (myrank == 0) clocks_gettime(&tic);
  struct engine e;
  engine_init(&e, &s, &params, nr_nodes, myrank, engine_policies, talking);
  if (talking) {
    clocks_gettime(&toc);
    message("engine_init took %.3f %s.", clocks_diff(&tic, &toc),
            clocks_getunit());
    fflush(stdout);
  }

#ifdef WITH_MPI
  /* Split the space. */
  engine_split(&e, &initial_partition);
  engine_redistribute(&e);
#endif

  int with_outputs = 1;
  if (with_outputs) {
    /* Write the state of the system before starting time integration. */
    if (myrank == 0) clocks_gettime(&tic);
#if defined(WITH_MPI)
#if defined(HAVE_PARALLEL_HDF5)
    write_output_parallel(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD,
                          MPI_INFO_NULL);
#else
    write_output_serial(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD,
                        MPI_INFO_NULL);
#endif
#else
    write_output_single(&e, &us);
#endif
    if (myrank == 0 && verbose) {
      clocks_gettime(&toc);
      message("writing particle properties took %.3f %s.",
              clocks_diff(&tic, &toc), clocks_getunit());
      fflush(stdout);
    }
  }

/* Init the runner history. */
#ifdef HIST
  for (k = 0; k < runner_hist_N; k++) runner_hist_bins[k] = 0;
#endif

  /* Get some info to the user. */
  if (myrank == 0) {
    message(
        "Running on %lld gas particles and %lld DM particles until t=%.3e with "
        "%i threads and %i queues (dt_min=%.3e, dt_max=%.3e)...",
        N_total[0], N_total[1], e.timeEnd, e.nr_threads, e.sched.nr_queues,
        e.dt_min, e.dt_max);
    fflush(stdout);
  }

  /* Initialise the particles */
  engine_init_particles(&e);

  /* Legend */
  if (myrank == 0)
    printf("# %6s %14s %14s %10s %10s %16s [%s]\n",
	   "Step",  "Time",  "Time-step",  "Updates", "g-Updates",
	   "Wall-clock time", clocks_getunit());

  /* Let loose a runner on the space. */
  for (int j = 0; !engine_is_done(&e); j++) {

/* Repartition the space amongst the nodes? */
#ifdef WITH_MPI
    if (j % 100 == 2) e.forcerepart = reparttype;
#endif

    timers_reset(timers_mask_all);
#ifdef COUNTER
    for (k = 0; k < runner_counter_count; k++) runner_counter[k] = 0;
#endif

    /* Take a step. */
    engine_step(&e);

    if (with_outputs && j % 100 == 0) {

      if (myrank == 0) clocks_gettime(&tic);
#if defined(WITH_MPI)
#if defined(HAVE_PARALLEL_HDF5)
      write_output_parallel(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD,
                            MPI_INFO_NULL);
#else
      write_output_serial(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD,
                          MPI_INFO_NULL);
#endif
#else
      write_output_single(&e, &us);
#endif
      if (myrank == 0 && verbose) {
        clocks_gettime(&toc);
        message("writing particle properties took %.3f %s.",
                clocks_diff(&tic, &toc), clocks_getunit());
        fflush(stdout);
      }
    }

    /* Dump the task data using the given frequency. */
    if (dump_tasks && (dump_tasks == 1 || j % dump_tasks == 1)) {
#ifdef WITH_MPI

      /* Make sure output file is empty, only on one rank. */
      char dumpfile[30];
      sprintf(dumpfile, "thread_info_MPI-step%d.dat", j);
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

          fprintf(file_thread, " %03i 0 0 0 0 %lli 0 0 0 0\n", myrank,
                  e.tic_step);
          int count = 0;
          for (int l = 0; l < e.sched.nr_tasks; l++)
            if (!e.sched.tasks[l].skip && !e.sched.tasks[l].implicit) {
              fprintf(file_thread, " %03i %i %i %i %i %lli %lli %i %i %i\n",
                      myrank, e.sched.tasks[l].rid, e.sched.tasks[l].type,
                      e.sched.tasks[l].subtype, (e.sched.tasks[l].cj == NULL),
                      e.sched.tasks[l].tic, e.sched.tasks[l].toc,
                      (e.sched.tasks[l].ci != NULL) ? e.sched.tasks[l].ci->count
                                                    : 0,
                      (e.sched.tasks[l].cj != NULL) ? e.sched.tasks[l].cj->count
                                                    : 0,
                      e.sched.tasks[l].flags);
              fflush(stdout);
              count++;
            }
          message("rank %d counted %d tasks", myrank, count);

          fclose(file_thread);
        }

        /* And we wait for all to synchronize. */
        MPI_Barrier(MPI_COMM_WORLD);
      }

#else
      char dumpfile[30];
      sprintf(dumpfile, "thread_info-step%d.dat", j);
      FILE *file_thread;
      file_thread = fopen(dumpfile, "w");
      for (int l = 0; l < e.sched.nr_tasks; l++)
        if (!e.sched.tasks[l].skip && !e.sched.tasks[l].implicit)
          fprintf(
              file_thread, " %i %i %i %i %lli %lli %i %i\n",
              e.sched.tasks[l].rid, e.sched.tasks[l].type,
              e.sched.tasks[l].subtype, (e.sched.tasks[l].cj == NULL),
              e.sched.tasks[l].tic, e.sched.tasks[l].toc,
              (e.sched.tasks[l].ci == NULL) ? 0 : e.sched.tasks[l].ci->count,
              (e.sched.tasks[l].cj == NULL) ? 0 : e.sched.tasks[l].cj->count);
      fclose(file_thread);
#endif
    }
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

  if (with_outputs) {

    if (myrank == 0) clocks_gettime(&tic);
/* Write final output. */
#if defined(WITH_MPI)
#if defined(HAVE_PARALLEL_HDF5)
    write_output_parallel(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD,
                          MPI_INFO_NULL);
#else
    write_output_serial(&e, &us, myrank, nr_nodes, MPI_COMM_WORLD,
                        MPI_INFO_NULL);
#endif
#else
    write_output_single(&e, &us);
#endif
    if (myrank == 0 && verbose) {
      clocks_gettime(&toc);
      message("writing particle properties took %.3f %s.",
              clocks_diff(&tic, &toc), clocks_getunit());
      fflush(stdout);
    }
  }

#ifdef WITH_MPI
  if (MPI_Finalize() != MPI_SUCCESS)
    error("call to MPI_Finalize failed with error %i.", res);
#endif

  /* Say goodbye. */
  if (myrank == 0) message("done. Bye.");

  /* All is calm, all is bright. */
  return 0;
}
