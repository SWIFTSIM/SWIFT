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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <fenv.h>

/* Conditional headers. */
#ifdef HAVE_LIBZ
#include <zlib.h>
#endif

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
 * @brief Main routine that loads a few particles and generates some output.
 *
 */
int main(int argc, char *argv[]) {

  int c, icount, periodic = 1;
  size_t Ngas = 0, Ngpart = 0;
  long long N_total[2] = {0, 0};
  int nr_threads = 1, nr_queues = -1;
  int dump_tasks = 0;
  int data[2];
  double dim[3] = {1.0, 1.0, 1.0}, shift[3] = {0.0, 0.0, 0.0};
  double h_max = -1.0, scaling = 1.0;
  double time_end = DBL_MAX;
  struct part *parts = NULL;
  struct gpart *gparts = NULL;
  struct space s;
  struct engine e;
  struct UnitSystem us;
  struct clocks_time tic, toc;
  char ICfileName[200] = "";
  char dumpfile[30];
  float dt_max = 0.0f, dt_min = 0.0f;
  int nr_nodes = 1, myrank = 0;
  FILE *file_thread;
  int with_outputs = 1;
  int with_external_gravity = 0;
  int with_self_gravity = 0;
  int engine_policies = 0;
  int verbose = 0, talking = 0;
  unsigned long long cpufreq = 0;

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

  /* Choke on FP-exceptions. */
  // feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );

  /* Initialize CPU frequency, this also starts time. */
  clocks_set_cpufreq(cpufreq);

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
  if (myrank == 0) message("MPI is up and running with %i node(s).", nr_nodes);
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
     * Otherwise,
     * we may be scheduled elsewhere between the two times.
     */
    cpu_set_t affinity;
    CPU_ZERO(&affinity);
    CPU_SET(sched_getcpu(), &affinity);
    if (sched_setaffinity(0, sizeof(cpu_set_t), &affinity) != 0) {
      message("failed to set entry thread's affinity");
    } else {
      message("set entry thread's affinity");
    }
  }
#endif

  /* Init the space. */
  bzero(&s, sizeof(struct space));

  /* Parse the options */
  while ((c = getopt(argc, argv, "a:c:d:e:f:gGh:m:oP:q:R:s:t:v:w:y:z:")) != -1)
    switch (c) {
      case 'a':
        if (sscanf(optarg, "%lf", &scaling) != 1)
          error("Error parsing cutoff scaling.");
        if (myrank == 0) message("scaling cutoff by %.3f.", scaling);
        fflush(stdout);
        break;
      case 'c':
        if (sscanf(optarg, "%lf", &time_end) != 1)
          error("Error parsing final time.");
        if (myrank == 0) message("time_end set to %.3e.", time_end);
        fflush(stdout);
        break;
      case 'd':
        if (sscanf(optarg, "%f", &dt_min) != 1)
          error("Error parsing minimal timestep.");
        if (myrank == 0) message("dt_min set to %e.", dt_min);
        fflush(stdout);
        break;
      case 'e':
        if (sscanf(optarg, "%f", &dt_max) != 1)
          error("Error parsing maximal timestep.");
        if (myrank == 0) message("dt_max set to %e.", dt_max);
        fflush(stdout);
        break;
      case 'f':
        if (!strcpy(ICfileName, optarg)) error("Error parsing IC file name.");
        break;
      case 'g':
        with_external_gravity = 1;
        break;
      case 'G':
        with_self_gravity = 1;
        break;
      case 'h':
        if (sscanf(optarg, "%llu", &cpufreq) != 1)
          error("Error parsing CPU frequency.");
        if (myrank == 0) message("CPU frequency set to %llu.", cpufreq);
        fflush(stdout);
        break;
      case 'm':
        if (sscanf(optarg, "%lf", &h_max) != 1) error("Error parsing h_max.");
        if (myrank == 0) message("maximum h set to %e.", h_max);
        fflush(stdout);
        break;
      case 'o':
        with_outputs = 0;
        break;
      case 'P':
/* Partition type is one of "g", "m", "w", or "v"; "g" can be
 * followed by three numbers defining the grid. */
#ifdef WITH_MPI
        switch (optarg[0]) {
          case 'g':
            initial_partition.type = INITPART_GRID;
            if (strlen(optarg) > 2) {
              if (sscanf(optarg, "g %i %i %i", &initial_partition.grid[0],
                         &initial_partition.grid[1],
                         &initial_partition.grid[2]) != 3)
                error("Error parsing grid.");
            }
            break;
#ifdef HAVE_METIS
          case 'm':
            initial_partition.type = INITPART_METIS_NOWEIGHT;
            break;
          case 'w':
            initial_partition.type = INITPART_METIS_WEIGHT;
            break;
#endif
          case 'v':
            initial_partition.type = INITPART_VECTORIZE;
            break;
        }
#endif
        break;
      case 'q':
        if (sscanf(optarg, "%d", &nr_queues) != 1)
          error("Error parsing number of queues.");
        break;
      case 'R':
/* Repartition type "n", "b", "v", "e" or "x".
 * Note only none is available without METIS. */
#ifdef WITH_MPI
        switch (optarg[0]) {
          case 'n':
            reparttype = REPART_NONE;
            break;
#ifdef HAVE_METIS
          case 'b':
            reparttype = REPART_METIS_BOTH;
            break;
          case 'v':
            reparttype = REPART_METIS_VERTEX;
            break;
          case 'e':
            reparttype = REPART_METIS_EDGE;
            break;
          case 'x':
            reparttype = REPART_METIS_VERTEX_EDGE;
            break;
#endif
        }
#endif
        break;
      case 's':
        if (sscanf(optarg, "%lf %lf %lf", &shift[0], &shift[1], &shift[2]) != 3)
          error("Error parsing shift.");
        if (myrank == 0)
          message("will shift parts by [ %.3f %.3f %.3f ].", shift[0], shift[1],
                  shift[2]);
        break;
      case 't':
        if (sscanf(optarg, "%d", &nr_threads) != 1)
          error("Error parsing number of threads.");
        break;
      case 'v':
        /* verbose = 1: MPI rank 0 writes
           verbose = 2: all MPI ranks write */
        if (sscanf(optarg, "%d", &verbose) != 1)
          error("Error parsing verbosity level.");
        break;
      case 'w':
        if (sscanf(optarg, "%d", &space_subsize) != 1)
          error("Error parsing sub size.");
        if (myrank == 0) message("sub size set to %i.", space_subsize);
        break;
      case 'y':
        if (sscanf(optarg, "%d", &dump_tasks) != 1)
          error("Error parsing dump_tasks (-y)");
        break;
      case 'z':
        if (sscanf(optarg, "%d", &space_splitsize) != 1)
          error("Error parsing split size.");
        if (myrank == 0) message("split size set to %i.", space_splitsize);
        break;
      case '?':
        error("Unknown option.");
        break;
    }

#ifdef WITH_MPI
  if (myrank == 0) {
    message("Running with %i thread(s) per node.", nr_threads);
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
  if (myrank == 0) message("Running with %i thread(s).", nr_threads);
#endif

  /* How large are the parts? */
  if (myrank == 0) {
    message("sizeof(struct part) is %li bytes.", (long int)sizeof(struct part));
    message("sizeof(struct xpart) is %li bytes.",
            (long int)sizeof(struct xpart));
    message("sizeof(struct gpart) is %li bytes.",
            (long int)sizeof(struct gpart));
  }

  /* Initialize unit system */
  initUnitSystem(&us);
  if (myrank == 0) {
    message("Unit system: U_M = %e g.", us.UnitMass_in_cgs);
    message("Unit system: U_L = %e cm.", us.UnitLength_in_cgs);
    message("Unit system: U_t = %e s.", us.UnitTime_in_cgs);
    message("Unit system: U_I = %e A.", us.UnitCurrent_in_cgs);
    message("Unit system: U_T = %e K.", us.UnitTemperature_in_cgs);
    message("Density units: %e a^%f h^%f.",
            conversionFactor(&us, UNIT_CONV_DENSITY),
            aFactor(&us, UNIT_CONV_DENSITY), hFactor(&us, UNIT_CONV_DENSITY));
    message("Entropy units: %e a^%f h^%f.",
            conversionFactor(&us, UNIT_CONV_ENTROPY),
            aFactor(&us, UNIT_CONV_ENTROPY), hFactor(&us, UNIT_CONV_ENTROPY));
  }

  /* Report CPU frequency. */
  if (myrank == 0) {
    cpufreq = clocks_get_cpufreq();
    message("CPU frequency used for tick conversion: %llu Hz", cpufreq);
  }

  /* Check whether an IC file has been provided */
  if (strcmp(ICfileName, "") == 0)
    error("An IC file name must be provided via the option -f");

  /* Read particles and space information from (GADGET) IC */

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
    message("reading particle properties took %.3f %s.",
            clocks_diff(&tic, &toc), clocks_getunit());
    fflush(stdout);
  }

#if defined(WITH_MPI)
  long long N_long[2] = {Ngas, Ngpart};
  MPI_Reduce(&N_long, &N_total, 2, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  N_total[1] -= N_total[0];
  if (myrank == 0)
    message("Read %lld gas particles and %lld DM particles from the ICs",
            N_total[0], N_total[1]);
#else
  N_total[0] = Ngas;
  N_total[1] = Ngpart - Ngas;
  message("Read %lld gas particles and %lld DM particles from the ICs",
          N_total[0], N_total[1]);
#endif

  /* MATTHIEU: Temporary fix to preserve master */
  if (!with_external_gravity && !with_self_gravity) {
    free(gparts);
    gparts = NULL;
    for (size_t k = 0; k < Ngas; ++k) parts[k].gpart = NULL;
    Ngpart = 0;
#if defined(WITH_MPI)
    N_long[0] = Ngas;
    N_long[1] = Ngpart;
    MPI_Reduce(&N_long, &N_total, 2, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0)
      message(
          "AFTER FIX: Read %lld gas particles and %lld DM particles from the "
          "ICs",
          N_total[0], N_total[1]);
#else
    N_total[0] = Ngas;
    N_total[1] = Ngpart;
    message(
        "AFTER FIX: Read %lld gas particles and %lld DM particles from the ICs",
        N_total[0], N_total[1]);
#endif
  }
  /* MATTHIEU: End temporary fix */

  /* Apply h scaling */
  if (scaling != 1.0)
    for (size_t k = 0; k < Ngas; k++) parts[k].h *= scaling;

  /* Apply shift */
  if (shift[0] != 0 || shift[1] != 0 || shift[2] != 0) {
    for (size_t k = 0; k < Ngas; k++) {
      parts[k].x[0] += shift[0];
      parts[k].x[1] += shift[1];
      parts[k].x[2] += shift[2];
    }
    for (size_t k = 0; k < Ngpart; k++) {
      gparts[k].x[0] += shift[0];
      gparts[k].x[1] += shift[1];
      gparts[k].x[2] += shift[2];
    }
  }

  /* Set default number of queues. */
  if (nr_queues < 0) nr_queues = nr_threads;

  /* How vocal are we ? */
  talking = (verbose == 1 && myrank == 0) || (verbose == 2);

  /* Initialize the space with this data. */
  if (myrank == 0) clocks_gettime(&tic);
  space_init(&s, dim, parts, gparts, Ngas, Ngpart, periodic, h_max,
             myrank == 0);
  if (myrank == 0 && verbose) {
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
    message("maximum depth is %d.", s.maxdepth);
    // message( "cutoffs in [ %g %g ]." , s.h_min , s.h_max ); fflush(stdout);
  }

  /* Verify that each particle is in it's proper cell. */
  if (myrank == 0) {
    icount = 0;
    space_map_cells_pre(&s, 0, &map_cellcheck, &icount);
    message("map_cellcheck picked up %i parts.", icount);
  }

  if (myrank == 0) {
    data[0] = s.maxdepth;
    data[1] = 0;
    space_map_cells_pre(&s, 0, &map_maxdepth, data);
    message("nr of cells at depth %i is %i.", data[0], data[1]);
  }

  /* Construct the engine policy */
  engine_policies = ENGINE_POLICY | engine_policy_steal | engine_policy_hydro;
  if (with_external_gravity) engine_policies |= engine_policy_external_gravity;
  if (with_self_gravity) engine_policies |= engine_policy_self_gravity;

  /* Initialize the engine with this space. */
  if (myrank == 0) clocks_gettime(&tic);
  if (myrank == 0) message("nr_nodes is %i.", nr_nodes);
  engine_init(&e, &s, dt_max, nr_threads, nr_queues, nr_nodes, myrank,
              engine_policies, 0, time_end, dt_min, dt_max, talking);
  if (myrank == 0 && verbose) {
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

  if (with_outputs) {
    /* Write the state of the system as it is before starting time integration.
     */
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

  if (myrank == 0) {
    message(
        "Running on %lld gas particles and %lld DM particles until t=%.3e with "
        "%i threads and %i queues (dt_min=%.3e, dt_max=%.3e)...",
        N_total[0], N_total[1], time_end, e.nr_threads, e.sched.nr_queues,
        e.dt_min, e.dt_max);
    fflush(stdout);
  }

  /* Initialise the particles */
  engine_init_particles(&e);

  /* Legend */
  if (myrank == 0)
    printf(
        "# Step  Time  time-step  Number of updates  Number of updates "
        "CPU Wall-clock time [%s]\n",
        clocks_getunit());

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
      sprintf(dumpfile, "thread_info_MPI-step%d.dat", j);
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
      sprintf(dumpfile, "thread_info-step%d.dat", j);
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

    /* Dump a line of aggregate output. */
    /*     if (myrank == 0) { */
    /*       printf("%i %e %.16e %.16e %.16e %.3e %.3e %i %.3e %.3e", j, e.time,
     */
    /*              e.ekin + e.epot, e.ekin, e.epot, e.dt, e.dt_step,
     * e.count_step, */
    /*              e.dt_min, e.dt_max); */
    /*       for (k = 0; k < timer_count; k++) */
    /*         printf(" %.3f", clocks_from_ticks(timers[k]); */
    /*       printf("\n"); */
    /*       fflush(stdout); */
    /*     } */

    /* if (myrank == 0) { */
    /*   printf("%i %e", j, e.time); */
    /*   printf(" %.3f", clocks_from_ticks(timers[timer_count - 1]); */
    /*   printf("\n"); */
    /*   fflush(stdout); */
    /* } */
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
