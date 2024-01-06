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

/* Config parameters. */
#include <config.h>

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <sched.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

/* MPI headers. */
#ifdef WITH_MPI

#include <mpi.h>
#endif

#ifdef HAVE_LIBNUMA
#include <numa.h>
#include <numaif.h>
#endif

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "active.h"
#include "atomic.h"
#include "black_holes_properties.h"
#include "cell.h"
#include "chemistry.h"
#include "clocks.h"
#include "cooling.h"
#include "cooling_properties.h"
#include "cosmology.h"
#include "csds.h"
#include "csds_io.h"
#include "cycle.h"
#include "debug.h"
#include "equation_of_state.h"
#include "error.h"
#include "extra_io.h"
#include "feedback.h"
#include "fof.h"
#include "forcing.h"
#include "gravity.h"
#include "gravity_cache.h"
#include "hydro.h"
#include "lightcone/lightcone.h"
#include "lightcone/lightcone_array.h"
#include "line_of_sight.h"
#include "map.h"
#include "memuse.h"
#include "minmax.h"
#include "mpiuse.h"
#include "multipole_struct.h"
#include "neutrino.h"
#include "neutrino_properties.h"
#include "output_list.h"
#include "output_options.h"
#include "partition.h"
#include "potential.h"
#include "power_spectrum.h"
#include "pressure_floor.h"
#include "profiler.h"
#include "proxy.h"
#include "restart.h"
#include "rt_properties.h"
#include "runner.h"
#include "sink_properties.h"
#include "sort_part.h"
#include "star_formation.h"
#include "star_formation_logger.h"
#include "stars_io.h"
#include "statistics.h"
#include "timers.h"
#include "tools.h"
#include "units.h"
#include "velociraptor_interface.h"

const char *engine_policy_names[] = {"none",
                                     "rand",
                                     "steal",
                                     "keep",
                                     "block",
                                     "cpu tight",
                                     "mpi",
                                     "numa affinity",
                                     "hydro",
                                     "self gravity",
                                     "external gravity",
                                     "cosmological integration",
                                     "drift everything",
                                     "reconstruct multi-poles",
                                     "temperature",
                                     "cooling",
                                     "stars",
                                     "structure finding",
                                     "star formation",
                                     "feedback",
                                     "black holes",
                                     "fof search",
                                     "time-step limiter",
                                     "time-step sync",
                                     "csds",
                                     "line of sight",
                                     "sink",
                                     "rt",
                                     "power spectra"};

const int engine_default_snapshot_subsample[swift_type_count] = {0};

/** The rank of the engine as a global variable (for messages). */
int engine_rank;

/** The current step of the engine as a global variable (for messages). */
int engine_current_step;

/**
 * @brief Link a density/force task to a cell.
 *
 * @param e The #engine.
 * @param l A pointer to the #link, will be modified atomically.
 * @param t The #task.
 *
 * @return The new #link pointer.
 */
void engine_addlink(struct engine *e, struct link **l, struct task *t) {

#ifdef SWIFT_DEBUG_CHECKS
  if (t == NULL) {
    error("Trying to link NULL task.");
  }
#endif

  /* Get the next free link. */
  const size_t ind = atomic_inc(&e->nr_links);
  if (ind >= e->size_links) {
    error(
        "Link table overflow. Increase the value of "
        "`Scheduler:links_per_tasks`.");
  }
  struct link *res = &e->links[ind];

  /* Set it atomically. */
  res->t = t;
  res->next = atomic_swap(l, res);
}

/**
 * @brief Repartition the cells amongst the nodes.
 *
 * @param e The #engine.
 */
void engine_repartition(struct engine *e) {

#if defined(WITH_MPI) && (defined(HAVE_PARMETIS) || defined(HAVE_METIS))

  ticks tic = getticks();

#ifdef SWIFT_DEBUG_CHECKS
  /* Be verbose about this. */
  if (e->nodeID == 0 || e->verbose) message("repartitioning space");
  fflush(stdout);

  /* Check that all cells have been drifted to the current time */
  space_check_drift_point(e->s, e->ti_current, /*check_multipoles=*/0);
#endif

  /* Clear the repartition flag. */
  e->forcerepart = 0;

  /* Nothing to do if only using a single node. Also avoids METIS
   * bug that doesn't handle this case well. */
  if (e->nr_nodes == 1) return;

  /* Generate the fixed costs include file. */
  if (e->step > 3 && e->reparttype->trigger <= 1.f) {
    task_dump_stats("partition_fixed_costs.h", e,
                    /* task_dump_threshold = */ 0.f,
                    /* header = */ 1, /* allranks = */ 1);
  }

  /* Do the repartitioning. */
  partition_repartition(e->reparttype, e->nodeID, e->nr_nodes, e->s,
                        e->sched.tasks, e->sched.nr_tasks);

  /* Partitioning requires copies of the particles, so we need to reduce the
   * memory in use to the minimum, we can free the sorting indices and the
   * tasks as these will be regenerated at the next rebuild. Also the foreign
   * particle arrays can go as these will be regenerated in proxy exchange. */

  /* Sorting indices. */
  if (e->s->cells_top != NULL) space_free_cells(e->s);

  /* Report the time spent in the different task categories */
  if (e->verbose) scheduler_report_task_times(&e->sched, e->nr_threads);

  /* Task arrays. */
  scheduler_free_tasks(&e->sched);

  /* Foreign parts. (no need to nullify the cell pointers as the cells
   * will be regenerated) */
  space_free_foreign_parts(e->s, /*clear_cell_pointers=*/0);

  /* Now comes the tricky part: Exchange particles between all nodes.
     This is done in two steps, first allreducing a matrix of
     how many particles go from where to where, then re-allocating
     the parts array, and emitting the sends and receives.
     Finally, the space, tasks, and proxies need to be rebuilt. */

  /* Redistribute the particles between the nodes. */
  engine_redistribute(e);

  /* Make the proxies. */
  engine_makeproxies(e);

  /* Tell the engine it should re-build whenever possible */
  e->forcerebuild = 1;

  /* Flag that a repartition has taken place */
  e->step_props |= engine_step_prop_repartition;

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  if (e->reparttype->type != REPART_NONE)
    error("SWIFT was not compiled with MPI and METIS or ParMETIS support.");

  /* Clear the repartition flag. */
  e->forcerepart = 0;
#endif
}

/**
 * @brief Decide whether trigger a repartition the cells amongst the nodes.
 *
 * @param e The #engine.
 */
void engine_repartition_trigger(struct engine *e) {

#ifdef WITH_MPI

  const ticks tic = getticks();
  static int opened = 0;
  if (e->restarting) opened = 1;

  /* Do nothing if there have not been enough steps since the last repartition
   * as we don't want to repeat this too often or immediately after a
   * repartition step, or also immediately on restart. We check all this
   * even when we are not repartitioning as the balance logs can still be
   * interesting. */
  if (e->step - e->last_repartition >= 2 && !e->restarting) {

    /* If we have fixed costs available and this is step 2 or we are forcing
     * repartitioning then we do a forced fixed costs repartition regardless. */
    int forced = 0;
    if (e->reparttype->type != REPART_NONE) {
      if (e->reparttype->trigger > 1 ||
          (e->step == 2 && e->reparttype->use_fixed_costs)) {
        if (e->reparttype->trigger > 1) {
          if ((e->step % (int)e->reparttype->trigger) == 0) e->forcerepart = 1;
        } else {
          e->forcerepart = 1;
        }
        e->reparttype->use_ticks = 0;
        forced = 1;
      }
    }

    /* We only check the CPU loads when we have processed a significant number
     * of all particles as we require all tasks to have timings or are
     * interested in the various balances logs. */
    if ((e->updates > 1 &&
         e->updates >= e->total_nr_parts * e->reparttype->minfrac) ||
        (e->g_updates > 1 &&
         e->g_updates >= e->total_nr_gparts * e->reparttype->minfrac)) {

      /* Are we using the task tick timings or fixed costs? */
      if (e->reparttype->use_fixed_costs > 1) {
        e->reparttype->use_ticks = 0;
      } else {
        e->reparttype->use_ticks = 1;
      }

      /* Get the resident size of the process for the memory logs. */
      long size, resident, shared, text, library, data, dirty;
      memuse_use(&size, &resident, &shared, &text, &data, &library, &dirty);

      /* Gather it together with the CPU times used by the tasks in the last
       * step (and the deadtime fraction, for output). */
      double timemem[4] = {
          e->usertime_last_step, e->systime_last_step, (double)resident,
          e->local_deadtime / (e->nr_threads * e->wallclock_time)};
      double timemems[e->nr_nodes * 4];
      MPI_Gather(&timemem, 4, MPI_DOUBLE, timemems, 4, MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
      if (e->nodeID == 0) {

        /* Get the range and mean of the two CPU times and memory. */
        double umintime = timemems[0];
        double umaxtime = timemems[0];

        double smintime = timemems[1];
        double smaxtime = timemems[1];

        double minmem = timemems[2];
        double maxmem = timemems[2];

        double tmintime = umintime + smintime;
        double tmaxtime = umaxtime + smaxtime;

        double usum = timemems[0];
        double ssum = timemems[1];
        double tsum = usum + ssum;

        double msum = timemems[2];

        for (int k = 4; k < e->nr_nodes * 4; k += 4) {
          if (timemems[k] > umaxtime) umaxtime = timemems[k];
          if (timemems[k] < umintime) umintime = timemems[k];

          if (timemems[k + 1] > smaxtime) smaxtime = timemems[k + 1];
          if (timemems[k + 1] < smintime) smintime = timemems[k + 1];

          double total = timemems[k] + timemems[k + 1];
          if (total > tmaxtime) tmaxtime = total;
          if (total < tmintime) tmintime = total;

          usum += timemems[k];
          ssum += timemems[k + 1];
          tsum += total;

          if (timemems[k + 2] > maxmem) maxmem = timemems[k + 2];
          if (timemems[k + 2] < minmem) minmem = timemems[k + 2];
          msum += timemems[k + 2];
        }
        double umean = usum / (double)e->nr_nodes;
        double smean = ssum / (double)e->nr_nodes;
        double tmean = tsum / (double)e->nr_nodes;
        double mmean = msum / (double)e->nr_nodes;

        /* Are we out of balance and need to repartition? */
        /* ---------------------------------------------- */
        double abs_trigger = fabs(e->reparttype->trigger);
        double balance = (umaxtime - umintime) / umean;
        if (e->reparttype->type != REPART_NONE) {

          /* When forced we don't care about the balance. */
          if (!forced) {
            if (balance > abs_trigger) {
              if (e->verbose)
                message("trigger fraction %.3f > %.3f will repartition",
                        balance, abs_trigger);
              e->forcerepart = 1;
            } else {
              if (e->verbose && e->reparttype->type != REPART_NONE)
                message("trigger fraction %.3f =< %.3f will not repartition",
                        balance, abs_trigger);
            }
          }

        } else {
          /* Not repartitioning, would that have been done otherwise? */
          if (e->verbose) {
            if (balance > abs_trigger) {
              message("trigger fraction %.3f > %.3f would have repartitioned",
                      balance, abs_trigger);
            } else {
              message(
                  "trigger fraction %.3f =< %.3f would not have repartitioned",
                  balance, abs_trigger);
            }
          }
        }

        /* Keep logs of all CPU times and resident memory size for debugging
         * load issues. */
        FILE *timelog = NULL;
        FILE *memlog = NULL;
        if (!opened) {
          timelog = fopen("rank_cpu_balance.log", "w");
          if (timelog == NULL)
            error("Could not create file 'rank_cpu_balance.log'.");
          fprintf(timelog, "# step rank user sys sum deadfrac\n");

          memlog = fopen("rank_memory_balance.log", "w");
          if (memlog == NULL)
            error("Could not create file 'rank_memory_balance.log'.");
          fprintf(memlog, "# step rank resident\n");

          opened = 1;
        } else {
          timelog = fopen("rank_cpu_balance.log", "a");
          if (timelog == NULL)
            error("Could not open file 'rank_cpu_balance.log' for writing.");
          memlog = fopen("rank_memory_balance.log", "a");
          if (memlog == NULL)
            error("Could not open file 'rank_memory_balance.log' for writing.");
        }

        for (int k = 0; k < e->nr_nodes * 4; k += 4) {

          fprintf(timelog, "%d %d %f %f %f %f\n", e->step, k / 4, timemems[k],
                  timemems[k + 1], timemems[k] + timemems[k + 1],
                  timemems[k + 3]);

          fprintf(memlog, "%d %d %ld\n", e->step, k / 4, (long)timemems[k + 2]);
        }

        fprintf(timelog, "# %d mean times: %f %f %f\n", e->step, umean, smean,
                tmean);
        if (abs_trigger > 1.f) abs_trigger = 0.f; /* Not relevant. */
        double systime = smean > 0. ? (smaxtime - smintime) / smean : 0.;
        double ttime = tmean > 0. ? (tmaxtime - tmintime) / tmean : 0.;
        fprintf(timelog,
                "# %d balance: %f, expected: %f (sys: %f, total: %f)\n",
                e->step, balance, abs_trigger, systime, ttime);

        fclose(timelog);

        fprintf(memlog, "# %d mean resident memory: %f, balance: %f\n", e->step,
                mmean, (maxmem - minmem) / mmean);
        fclose(memlog);
      }

      /* All nodes do this together, so send to other ranks. */
      MPI_Bcast(&e->forcerepart, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    /* Remember we did this. */
    if (e->forcerepart) e->last_repartition = e->step;
  }

  if (e->verbose)
    message("took %.3f %s", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#endif
}

/**
 * @brief Exchange cell structures with other nodes.
 *
 * @param e The #engine.
 */
void engine_exchange_cells(struct engine *e) {

#ifdef WITH_MPI

  const int with_gravity = e->policy & engine_policy_self_gravity;
  const ticks tic = getticks();

  /* Exchange the cell structure with neighbouring ranks. */
  proxy_cells_exchange(e->proxies, e->nr_proxies, e->s, with_gravity);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Exchanges the top-level multipoles between all the nodes
 * such that every node has a multipole for each top-level cell.
 *
 * @param e The #engine.
 */
void engine_exchange_top_multipoles(struct engine *e) {

#ifdef WITH_MPI

  ticks tic = getticks();

#ifdef SWIFT_DEBUG_CHECKS
  for (int i = 0; i < e->s->nr_cells; ++i) {
    const struct gravity_tensors *m = &e->s->multipoles_top[i];
    if (e->s->cells_top[i].nodeID == engine_rank) {
      if (m->m_pole.M_000 > 0.) {
        if (m->CoM[0] < 0. || m->CoM[0] > e->s->dim[0])
          error("Invalid multipole position in X");
        if (m->CoM[1] < 0. || m->CoM[1] > e->s->dim[1])
          error("Invalid multipole position in Y");
        if (m->CoM[2] < 0. || m->CoM[2] > e->s->dim[2])
          error("Invalid multipole position in Z");
      }
    } else {
      if (m->m_pole.M_000 != 0.) error("Non-zero mass for foreign m-pole");
      if (m->CoM[0] != 0.) error("Non-zero position in X for foreign m-pole");
      if (m->CoM[1] != 0.) error("Non-zero position in Y for foreign m-pole");
      if (m->CoM[2] != 0.) error("Non-zero position in Z for foreign m-pole");
      if (m->m_pole.num_gpart != 0)
        error("Non-zero gpart count in foreign m-pole");
    }
  }
#endif

  /* Each node (space) has constructed its own top-level multipoles.
   * We now need to make sure every other node has a copy of everything.
   *
   * We use our home-made reduction operation that simply performs a XOR
   * operation on the multipoles. Since only local multipoles are non-zero and
   * each multipole is only present once, the bit-by-bit XOR will
   * create the desired result.
   */
  int err = MPI_Allreduce(MPI_IN_PLACE, e->s->multipoles_top, e->s->nr_cells,
                          multipole_mpi_type, multipole_mpi_reduce_op,
                          MPI_COMM_WORLD);
  if (err != MPI_SUCCESS)
    mpi_error(err, "Failed to all-reduce the top-level multipoles.");

#ifdef SWIFT_DEBUG_CHECKS
  long long counter = 0;

  /* Let's check that what we received makes sense */
  for (int i = 0; i < e->s->nr_cells; ++i) {
    const struct gravity_tensors *m = &e->s->multipoles_top[i];
    counter += m->m_pole.num_gpart;
    if (m->m_pole.num_gpart < 0) {
      error("m->m_pole.num_gpart is negative: %lld", m->m_pole.num_gpart);
    }
    if (m->m_pole.M_000 > 0.) {
      if (m->CoM[0] < 0. || m->CoM[0] > e->s->dim[0])
        error("Invalid multipole position in X");
      if (m->CoM[1] < 0. || m->CoM[1] > e->s->dim[1])
        error("Invalid multipole position in Y");
      if (m->CoM[2] < 0. || m->CoM[2] > e->s->dim[2])
        error("Invalid multipole position in Z");
    }
  }
  if (counter != e->total_nr_gparts)
    error(
        "Total particles in multipoles inconsistent with engine.\n "
        "  counter = %lld, nr_gparts = %lld",
        counter, e->total_nr_gparts);
#endif

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

void engine_exchange_proxy_multipoles(struct engine *e) {

#ifdef WITH_MPI

  const ticks tic = getticks();

  /* Start by counting the number of cells to send and receive */
  int count_send_cells = 0;
  int count_recv_cells = 0;
  int count_send_requests = 0;
  int count_recv_requests = 0;

  /* Loop over the proxies. */
  for (int pid = 0; pid < e->nr_proxies; pid++) {

    /* Get a handle on the proxy. */
    const struct proxy *p = &e->proxies[pid];

    /* Now collect the number of requests associated */
    count_recv_requests += p->nr_cells_in;
    count_send_requests += p->nr_cells_out;

    /* And the actual number of things we are going to ship */
    for (int k = 0; k < p->nr_cells_in; k++)
      count_recv_cells += p->cells_in[k]->mpi.pcell_size;

    for (int k = 0; k < p->nr_cells_out; k++)
      count_send_cells += p->cells_out[k]->mpi.pcell_size;
  }

  /* Allocate the buffers for the packed data */
  struct gravity_tensors *buffer_send = NULL;
  if (swift_memalign("send_gravity_tensors", (void **)&buffer_send,
                     SWIFT_CACHE_ALIGNMENT,
                     count_send_cells * sizeof(struct gravity_tensors)) != 0)
    error("Unable to allocate memory for multipole transactions");

  struct gravity_tensors *buffer_recv = NULL;
  if (swift_memalign("recv_gravity_tensors", (void **)&buffer_recv,
                     SWIFT_CACHE_ALIGNMENT,
                     count_recv_cells * sizeof(struct gravity_tensors)) != 0)
    error("Unable to allocate memory for multipole transactions");

  /* Also allocate the MPI requests */
  const int count_requests = count_send_requests + count_recv_requests;
  MPI_Request *requests =
      (MPI_Request *)malloc(sizeof(MPI_Request) * count_requests);
  if (requests == NULL) error("Unable to allocate memory for MPI requests");

  int this_request = 0;
  int this_recv = 0;
  int this_send = 0;

  /* Loop over the proxies to issue the receives. */
  for (int pid = 0; pid < e->nr_proxies; pid++) {

    /* Get a handle on the proxy. */
    const struct proxy *p = &e->proxies[pid];

    for (int k = 0; k < p->nr_cells_in; k++) {

      const int num_elements = p->cells_in[k]->mpi.pcell_size;

      /* Receive everything */
      MPI_Irecv(&buffer_recv[this_recv], num_elements, multipole_mpi_type,
                p->cells_in[k]->nodeID, p->cells_in[k]->mpi.tag, MPI_COMM_WORLD,
                &requests[this_request]);

      /* Move to the next slot in the buffers */
      this_recv += num_elements;
      this_request++;
    }

    /* Loop over the proxies to issue the sends. */
    for (int k = 0; k < p->nr_cells_out; k++) {

      /* Number of multipoles in this cell hierarchy */
      const int num_elements = p->cells_out[k]->mpi.pcell_size;

      /* Let's pack everything recursively */
      cell_pack_multipoles(p->cells_out[k], &buffer_send[this_send]);

      /* Send everything (note the use of cells_in[0] to get the correct node
       * ID. */
      MPI_Isend(&buffer_send[this_send], num_elements, multipole_mpi_type,
                p->cells_in[0]->nodeID, p->cells_out[k]->mpi.tag,
                MPI_COMM_WORLD, &requests[this_request]);

      /* Move to the next slot in the buffers */
      this_send += num_elements;
      this_request++;
    }
  }

  /* Wait for all the requests to arrive home */
  MPI_Status *stats = (MPI_Status *)malloc(count_requests * sizeof(MPI_Status));
  int res;
  if ((res = MPI_Waitall(count_requests, requests, stats)) != MPI_SUCCESS) {
    for (int k = 0; k < count_requests; ++k) {
      char buff[MPI_MAX_ERROR_STRING];
      MPI_Error_string(stats[k].MPI_ERROR, buff, &res);
      message("request from source %i, tag %i has error '%s'.",
              stats[k].MPI_SOURCE, stats[k].MPI_TAG, buff);
    }
    error("Failed during waitall for multipole data.");
  }

  /* Let's now unpack the multipoles at the right place */
  this_recv = 0;
  for (int pid = 0; pid < e->nr_proxies; pid++) {

    /* Get a handle on the proxy. */
    const struct proxy *p = &e->proxies[pid];

    for (int k = 0; k < p->nr_cells_in; k++) {

      const int num_elements = p->cells_in[k]->mpi.pcell_size;

#ifdef SWIFT_DEBUG_CHECKS

      /* Check that the first element (top-level cell's multipole) matches what
       * we received */
      if (p->cells_in[k]->grav.multipole->m_pole.num_gpart !=
          buffer_recv[this_recv].m_pole.num_gpart)
        error("Current: M_000=%e num_gpart=%lld\n New: M_000=%e num_gpart=%lld",
              p->cells_in[k]->grav.multipole->m_pole.M_000,
              p->cells_in[k]->grav.multipole->m_pole.num_gpart,
              buffer_recv[this_recv].m_pole.M_000,
              buffer_recv[this_recv].m_pole.num_gpart);
#endif

      /* Unpack recursively */
      cell_unpack_multipoles(p->cells_in[k], &buffer_recv[this_recv]);

      /* Move to the next slot in the buffers */
      this_recv += num_elements;
    }
  }

  /* Free everything */
  free(stats);
  free(buffer_send);
  free(buffer_recv);
  free(requests);

  /* How much time did this take? */
  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Allocate memory for the foreign particles.
 *
 * We look into the proxies for cells that have tasks and count
 * the number of particles in these cells. We then allocate
 * memory and link all the cells that have tasks and all cells
 * deeper in the tree.
 *
 * When running FOF, we only need #gpart arrays so we restrict
 * the allocations to this particle type only
 *
 * @param e The #engine.
 * @param fof Are we allocating buffers just for FOF?
 */
void engine_allocate_foreign_particles(struct engine *e, const int fof) {

#ifdef WITH_MPI

  const int nr_proxies = e->nr_proxies;
  const int with_hydro = e->policy & engine_policy_hydro;
  const int with_stars = e->policy & engine_policy_stars;
  const int with_black_holes = e->policy & engine_policy_black_holes;
  struct space *s = e->s;
  ticks tic = getticks();

  /* Count the number of particles we need to import and re-allocate
     the buffer if needed. */
  size_t count_parts_in = 0, count_gparts_in = 0, count_sparts_in = 0,
         count_bparts_in = 0;
  for (int k = 0; k < nr_proxies; k++) {
    for (int j = 0; j < e->proxies[k].nr_cells_in; j++) {

      if (e->proxies[k].cells_in_type[j] & proxy_cell_type_hydro) {
        count_parts_in += cell_count_parts_for_tasks(e->proxies[k].cells_in[j]);
      }

      if (e->proxies[k].cells_in_type[j] & proxy_cell_type_gravity) {
        count_gparts_in +=
            cell_count_gparts_for_tasks(e->proxies[k].cells_in[j]);
      }

      /* For stars, we just use the numbers in the top-level cells */
      count_sparts_in +=
          e->proxies[k].cells_in[j]->stars.count + space_extra_sparts;

      /* For black holes, we just use the numbers in the top-level cells */
      count_bparts_in += e->proxies[k].cells_in[j]->black_holes.count;
    }
  }

  if (!with_hydro && count_parts_in)
    error(
        "Not running with hydro but about to receive gas particles in "
        "proxies!");

  if (e->verbose)
    message("Counting number of foreign particles took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  tic = getticks();

  /* Allocate space for the foreign particles we will receive */
  size_t old_size_parts_foreign = s->size_parts_foreign;
  if (!fof && count_parts_in > s->size_parts_foreign) {
    if (s->parts_foreign != NULL) swift_free("parts_foreign", s->parts_foreign);
    s->size_parts_foreign = engine_foreign_alloc_margin * count_parts_in;
    if (swift_memalign("parts_foreign", (void **)&s->parts_foreign, part_align,
                       sizeof(struct part) * s->size_parts_foreign) != 0)
      error("Failed to allocate foreign part data.");
  }

  /* Allocate space for the foreign particles we will receive */
  size_t old_size_gparts_foreign = s->size_gparts_foreign;
  if (count_gparts_in > s->size_gparts_foreign) {
    if (s->gparts_foreign != NULL)
      swift_free("gparts_foreign", s->gparts_foreign);
    s->size_gparts_foreign = engine_foreign_alloc_margin * count_gparts_in;
    if (swift_memalign("gparts_foreign", (void **)&s->gparts_foreign,
                       gpart_align,
                       sizeof(struct gpart) * s->size_gparts_foreign) != 0)
      error("Failed to allocate foreign gpart data.");
  }

  /* Allocate space for the foreign particles we will receive */
  size_t old_size_sparts_foreign = s->size_sparts_foreign;
  if (!fof && count_sparts_in > s->size_sparts_foreign) {
    if (s->sparts_foreign != NULL)
      swift_free("sparts_foreign", s->sparts_foreign);
    s->size_sparts_foreign = engine_foreign_alloc_margin * count_sparts_in;
    if (swift_memalign("sparts_foreign", (void **)&s->sparts_foreign,
                       spart_align,
                       sizeof(struct spart) * s->size_sparts_foreign) != 0)
      error("Failed to allocate foreign spart data.");
    bzero(s->sparts_foreign, s->size_sparts_foreign * sizeof(struct spart));
    for (size_t i = 0; i < s->size_sparts_foreign; ++i) {
      s->sparts_foreign[i].time_bin = time_bin_not_created;
      s->sparts_foreign[i].id = -43;
    }
  }

  /* Allocate space for the foreign particles we will receive */
  size_t old_size_bparts_foreign = s->size_bparts_foreign;
  if (!fof && count_bparts_in > s->size_bparts_foreign) {
    if (s->bparts_foreign != NULL)
      swift_free("bparts_foreign", s->bparts_foreign);
    s->size_bparts_foreign = engine_foreign_alloc_margin * count_bparts_in;
    if (swift_memalign("bparts_foreign", (void **)&s->bparts_foreign,
                       bpart_align,
                       sizeof(struct bpart) * s->size_bparts_foreign) != 0)
      error("Failed to allocate foreign bpart data.");
  }

  if (e->verbose) {
    message(
        "Allocating %zd/%zd/%zd/%zd foreign part/gpart/spart/bpart "
        "(%zd/%zd/%zd/%zd MB)",
        s->size_parts_foreign, s->size_gparts_foreign, s->size_sparts_foreign,
        s->size_bparts_foreign,
        s->size_parts_foreign * sizeof(struct part) / (1024 * 1024),
        s->size_gparts_foreign * sizeof(struct gpart) / (1024 * 1024),
        s->size_sparts_foreign * sizeof(struct spart) / (1024 * 1024),
        s->size_bparts_foreign * sizeof(struct bpart) / (1024 * 1024));

    if ((s->size_parts_foreign - old_size_parts_foreign) > 0 ||
        (s->size_gparts_foreign - old_size_gparts_foreign) > 0 ||
        (s->size_sparts_foreign - old_size_sparts_foreign) > 0 ||
        (s->size_bparts_foreign - old_size_bparts_foreign) > 0) {
      message(
          "Re-allocations %zd/%zd/%zd/%zd part/gpart/spart/bpart "
          "(%zd/%zd/%zd/%zd MB)",
          (s->size_parts_foreign - old_size_parts_foreign),
          (s->size_gparts_foreign - old_size_gparts_foreign),
          (s->size_sparts_foreign - old_size_sparts_foreign),
          (s->size_bparts_foreign - old_size_bparts_foreign),
          (s->size_parts_foreign - old_size_parts_foreign) *
              sizeof(struct part) / (1024 * 1024),
          (s->size_gparts_foreign - old_size_gparts_foreign) *
              sizeof(struct gpart) / (1024 * 1024),
          (s->size_sparts_foreign - old_size_sparts_foreign) *
              sizeof(struct spart) / (1024 * 1024),
          (s->size_bparts_foreign - old_size_bparts_foreign) *
              sizeof(struct bpart) / (1024 * 1024));
    }
  }

  /* Unpack the cells and link to the particle data. */
  struct part *parts = s->parts_foreign;
  struct gpart *gparts = s->gparts_foreign;
  struct spart *sparts = s->sparts_foreign;
  struct bpart *bparts = s->bparts_foreign;
  for (int k = 0; k < nr_proxies; k++) {
    for (int j = 0; j < e->proxies[k].nr_cells_in; j++) {

      if (!fof && e->proxies[k].cells_in_type[j] & proxy_cell_type_hydro) {

        const size_t count_parts =
            cell_link_foreign_parts(e->proxies[k].cells_in[j], parts);
        parts = &parts[count_parts];
      }

      if (e->proxies[k].cells_in_type[j] & proxy_cell_type_gravity) {

        const size_t count_gparts =
            cell_link_foreign_gparts(e->proxies[k].cells_in[j], gparts);
        gparts = &gparts[count_gparts];
      }

      if (!fof && with_stars) {

        /* For stars, we just use the numbers in the top-level cells */
        cell_link_sparts(e->proxies[k].cells_in[j], sparts);
        sparts = &sparts[e->proxies[k].cells_in[j]->stars.count +
                         space_extra_sparts];
      }

      if (!fof && with_black_holes) {

        /* For black holes, we just use the numbers in the top-level cells */
        cell_link_bparts(e->proxies[k].cells_in[j], bparts);
        bparts = &bparts[e->proxies[k].cells_in[j]->black_holes.count];
      }
    }
  }

  /* Update the counters */
  s->nr_parts_foreign = parts - s->parts_foreign;
  s->nr_gparts_foreign = gparts - s->gparts_foreign;
  s->nr_sparts_foreign = sparts - s->sparts_foreign;
  s->nr_bparts_foreign = bparts - s->bparts_foreign;

  if (e->verbose)
    message("Recursively linking foreign arrays took %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

void engine_do_tasks_count_mapper(void *map_data, int num_elements,
                                  void *extra_data) {

  const struct task *tasks = (struct task *)map_data;
  int *const global_counts = (int *)extra_data;

  /* Local accumulator copy */
  int local_counts[task_type_count + 1];
  for (int k = 0; k <= task_type_count; k++) local_counts[k] = 0;

  /* Add task counts locally */
  for (int k = 0; k < num_elements; k++) {
    if (tasks[k].skip)
      local_counts[task_type_count] += 1;
    else
      local_counts[(int)tasks[k].type] += 1;
  }

  /* Update the global counts */
  for (int k = 0; k <= task_type_count; k++) {
    if (local_counts[k]) atomic_add(global_counts + k, local_counts[k]);
  }
}

/**
 * @brief Prints the number of tasks in the engine
 *
 * @param e The #engine.
 */
void engine_print_task_counts(const struct engine *e) {

  const ticks tic = getticks();
  const struct scheduler *sched = &e->sched;
  const int nr_tasks = sched->nr_tasks;
  const struct task *const tasks = sched->tasks;

  /* Global tasks and cells when using MPI. */
#ifdef WITH_MPI
  if (e->nodeID == 0 && e->total_nr_tasks > 0)
    printf(
        "[%04i] %s engine_print_task_counts: System total: %lld,"
        " no. cells: %lld\n",
        e->nodeID, clocks_get_timesincestart(), e->total_nr_tasks,
        e->total_nr_cells);
  fflush(stdout);
#endif

  /* Report value that can be used to estimate the task_per_cells parameter. */
  float tasks_per_cell = (float)nr_tasks / (float)e->s->tot_cells;

#ifdef WITH_MPI
  message("Total = %d (per cell = %.2f)", nr_tasks, tasks_per_cell);

  /* And the system maximum on rank 0, only after first step, increase by our
   * margin to allow for some variation in repartitioning. */
  if (e->nodeID == 0 && e->total_nr_tasks > 0) {
    message("Total = %d (maximum per cell = %.2f)", nr_tasks,
            e->tasks_per_cell_max * engine_tasks_per_cell_margin);
  }

#else
  message("Total = %d (per cell = %.2f)", nr_tasks, tasks_per_cell);
#endif
  fflush(stdout);

  /* Count and print the number of each task type. */
  int counts[task_type_count + 1];
  for (int k = 0; k <= task_type_count; k++) counts[k] = 0;
  threadpool_map((struct threadpool *)&e->threadpool,
                 engine_do_tasks_count_mapper, (void *)tasks, nr_tasks,
                 sizeof(struct task), threadpool_auto_chunk_size, counts);

#ifdef WITH_MPI
  printf("[%04i] %s engine_print_task_counts: task counts are [ %s=%i",
         e->nodeID, clocks_get_timesincestart(), taskID_names[0], counts[0]);
#else
  printf("%s engine_print_task_counts: task counts are [ %s=%i",
         clocks_get_timesincestart(), taskID_names[0], counts[0]);
#endif
  for (int k = 1; k < task_type_count; k++)
    printf(" %s=%i", taskID_names[k], counts[k]);
  printf(" skipped=%i ]\n", counts[task_type_count]);
  fflush(stdout);
  message("nr_parts = %zu.", e->s->nr_parts);
  message("nr_gparts = %zu.", e->s->nr_gparts);
  message("nr_sink = %zu.", e->s->nr_sinks);
  message("nr_sparts = %zu.", e->s->nr_sparts);
  message("nr_bparts = %zu.", e->s->nr_bparts);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief if necessary, estimate the number of tasks required given
 *        the current tasks in use and the numbers of cells.
 *
 * If e->tasks_per_cell is set greater than 0.0 then that value is used
 * as the estimate of the average number of tasks per cell,
 * otherwise we attempt an estimate.
 *
 * @param e the #engine
 *
 * @return the estimated total number of tasks
 */
int engine_estimate_nr_tasks(const struct engine *e) {

  float tasks_per_cell = e->tasks_per_cell;
  if (tasks_per_cell > 0.0f) {
    if (e->verbose)
      message("tasks per cell given as: %.2f, so maximum tasks: %d",
              e->tasks_per_cell, (int)(e->s->tot_cells * tasks_per_cell));
    return (int)(e->s->tot_cells * tasks_per_cell);
  }

  /* Our guess differs depending on the types of tasks we are using, but we
   * basically use a formula <n1>*ntopcells + <n2>*(totcells - ntopcells).
   * Where <n1> is the expected maximum tasks per top-level/super cell, and
   * <n2> the expected maximum tasks for all other cells. These should give
   * a safe upper limit. */
  int n1 = 0;
  int n2 = 0;
  if (e->policy & engine_policy_hydro) {
    /* 2 self (density, force), 1 sort, 26/2 density pairs
       26/2 force pairs, 1 drift, 3 ghosts, 2 kicks, 1 time-step,
       1 end_force, 1 collect, 2 extra space
     */
    n1 += 38;
    n2 += 2;
#ifdef WITH_MPI
    n1 += 6;
#endif

#ifdef EXTRA_HYDRO_LOOP
    n1 += 15;
#ifdef WITH_MPI
    n1 += 2;
#endif
#endif
  }
  if (e->policy & engine_policy_timestep_limiter) {
    n1 += 18;
    n2 += 1;
  }
  if (e->policy & engine_policy_self_gravity) {
    n1 += 125;
    n2 += 8;
#ifdef WITH_MPI
    n2 += 2;
#endif
  }
  if (e->policy & engine_policy_external_gravity) {
    n1 += 2;
  }
  if (e->policy & engine_policy_cosmology) {
    n1 += 2;
  }
  if (e->policy & engine_policy_cooling) {
    /* Cooling task + extra space */
    n1 += 2;
  }
  if (e->policy & engine_policy_star_formation) {
    n1 += 1;
  }
  if (e->policy & engine_policy_stars) {
    /* 2 self (density, feedback), 1 sort, 26/2 density pairs
       26/2 feedback pairs, 1 drift, 3 ghosts, 2 kicks, 1 time-step,
       1 end_force, 2 extra space
     */
    n1 += 37;
    n2 += 2;
#ifdef WITH_MPI
    n1 += 6;
#endif
  }
  if (e->policy & engine_policy_sinks) {
    /* 1 drift, 2 kicks, 1 time-step, 1 sink formation */
    n1 += 5;
    if (e->policy & engine_policy_stars) {
      /* 1 star formation */
      n1 += 1;
    }
  }
  if (e->policy & engine_policy_fof) {
    n1 += 2;
  }
#if defined(WITH_CSDS)
  /* each cell logs its particles */
  if (e->policy & engine_policy_csds) {
    n1 += 1;
  }
#endif
  if (e->policy & engine_policy_rt) {
    /* gradient: 1 self + 13 pairs                   |   14
     * transport: 1 self + 13 pairs                  | + 14
     * implicits: in + out, transport_out            | +  3
     * others: ghost1, ghost2, tchem, cell advance   | +  4
     * sort, collect_times                           | +  2
     * 2 extra space                                 | +  2 */
    n1 += 39;
#ifdef WITH_MPI
    n1 += 2;
#endif
  }

#ifdef WITH_MPI

  /* We need fewer tasks per rank when using MPI, but we could have
   * imbalances, so we need to work using the locally active cells, not just
   * some equipartition amongst the nodes. Don't want to recurse the whole
   * cell tree, so just make a guess of the maximum possible total cells. */
  int ntop = 0;
  int ncells = 0;
  for (int k = 0; k < e->s->nr_cells; k++) {
    struct cell *c = &e->s->cells_top[k];

    /* Any cells with particles will have tasks (local & foreign). */
    int nparts = c->hydro.count + c->grav.count + c->stars.count;
    if (nparts > 0) {
      ntop++;
      ncells++;

      /* Count cell depth until we get below the parts per cell threshold. */
      int depth = 0;
      while (nparts > space_splitsize) {
        depth++;
        nparts /= 8;
        ncells += (1 << (depth * 3));
      }
    }
  }

  /* If no local cells, we are probably still initialising, so just keep
   * room for the top-level. */
  if (ncells == 0) {
    ntop = e->s->nr_cells;
    ncells = ntop;
  }
#else
  int ntop = e->s->nr_cells;
  int ncells = e->s->tot_cells;
#endif

  float ntasks = n1 * ntop + n2 * (ncells - ntop);
  if (ncells > 0) tasks_per_cell = ceilf(ntasks / ncells);

  if (tasks_per_cell < 1.0f) tasks_per_cell = 1.0f;
  if (e->verbose)
    message("tasks per cell estimated as: %.2f, maximum tasks: %d",
            tasks_per_cell, (int)(ncells * tasks_per_cell));

  return (int)(ncells * tasks_per_cell);
}

/**
 * @brief Rebuild the space and tasks.
 *
 * @param e The #engine.
 * @param repartitioned Did we just redistribute?
 * @param clean_smoothing_length_values Are we cleaning up the values of
 * the smoothing lengths before building the tasks ?
 */
void engine_rebuild(struct engine *e, const int repartitioned,
                    const int clean_smoothing_length_values) {

  const ticks tic = getticks();

  /* Clear the forcerebuild flag, whatever it was. */
  e->forcerebuild = 0;
  e->restarting = 0;

  /* Report the time spent in the different task categories */
  if (e->verbose && !repartitioned)
    scheduler_report_task_times(&e->sched, e->nr_threads);

  /* Give some breathing space */
  scheduler_free_tasks(&e->sched);

  /* Free the foreign particles to get more breathing space. */
#ifdef WITH_MPI
  if (e->free_foreign_when_rebuilding)
    space_free_foreign_parts(e->s, /*clear_cell_pointers=*/1);
#endif

  /* Re-build the space. */
  space_rebuild(e->s, repartitioned, e->verbose);

  /* Report the number of cells and memory */
  if (e->verbose)
    message(
        "Nr. of top-level cells: %d Nr. of local cells: %d memory use: %zd MB.",
        e->s->nr_cells, e->s->tot_cells,
        (e->s->nr_cells + e->s->tot_cells) * sizeof(struct cell) /
            (1024 * 1024));

  /* Report the number of multipoles and memory */
  if (e->verbose && (e->policy & engine_policy_self_gravity))
    message(
        "Nr. of top-level mpoles: %d Nr. of local mpoles: %d memory use: %zd "
        "MB.",
        e->s->nr_cells, e->s->tot_cells,
        (e->s->nr_cells + e->s->tot_cells) * sizeof(struct gravity_tensors) /
            (1024 * 1024));

  /* Report the number of particles and memory */
  if (e->verbose)
    message(
        "Space has memory for %zd/%zd/%zd/%zd/%zd part/gpart/spart/sink/bpart "
        "(%zd/%zd/%zd/%zd/%zd MB)",
        e->s->size_parts, e->s->size_gparts, e->s->size_sparts,
        e->s->size_sinks, e->s->size_bparts,
        e->s->size_parts * sizeof(struct part) / (1024 * 1024),
        e->s->size_gparts * sizeof(struct gpart) / (1024 * 1024),
        e->s->size_sparts * sizeof(struct spart) / (1024 * 1024),
        e->s->size_sinks * sizeof(struct sink) / (1024 * 1024),
        e->s->size_bparts * sizeof(struct bpart) / (1024 * 1024));

  if (e->verbose)
    message(
        "Space holds %zd/%zd/%zd/%zd/%zd part/gpart/spart/sink/bpart (fracs: "
        "%f/%f/%f/%f/%f)",
        e->s->nr_parts, e->s->nr_gparts, e->s->nr_sparts, e->s->nr_sinks,
        e->s->nr_bparts,
        e->s->nr_parts ? e->s->nr_parts / ((double)e->s->size_parts) : 0.,
        e->s->nr_gparts ? e->s->nr_gparts / ((double)e->s->size_gparts) : 0.,
        e->s->nr_sparts ? e->s->nr_sparts / ((double)e->s->size_sparts) : 0.,
        e->s->nr_sinks ? e->s->nr_sinks / ((double)e->s->size_sinks) : 0.,
        e->s->nr_bparts ? e->s->nr_bparts / ((double)e->s->size_bparts) : 0.);

  const ticks tic2 = getticks();

  /* Update the global counters of particles */
  long long num_particles[5] = {
      (long long)(e->s->nr_parts - e->s->nr_extra_parts),
      (long long)(e->s->nr_gparts - e->s->nr_extra_gparts),
      (long long)(e->s->nr_sparts - e->s->nr_extra_sparts),
      (long long)(e->s->nr_sinks - e->s->nr_extra_sinks),
      (long long)(e->s->nr_bparts - e->s->nr_extra_bparts)};
#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, num_particles, 5, MPI_LONG_LONG, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  e->total_nr_parts = num_particles[0];
  e->total_nr_gparts = num_particles[1];
  e->total_nr_sparts = num_particles[2];
  e->total_nr_sinks = num_particles[3];
  e->total_nr_bparts = num_particles[4];

#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &e->s->min_a_grav, 1, MPI_FLOAT, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &e->s->max_softening, 1, MPI_FLOAT, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &e->s->max_mpole_power,
                SELF_GRAVITY_MULTIPOLE_ORDER + 1, MPI_FLOAT, MPI_MAX,
                MPI_COMM_WORLD);
#endif

  /* Flag that there are no inhibited particles */
  e->nr_inhibited_parts = 0;
  e->nr_inhibited_gparts = 0;
  e->nr_inhibited_sparts = 0;
  e->nr_inhibited_sinks = 0;
  e->nr_inhibited_bparts = 0;

  if (e->verbose)
    message("updating particle counts took %.3f %s.",
            clocks_from_ticks(getticks() - tic2), clocks_getunit());

#ifdef SWIFT_DEBUG_CHECKS
  part_verify_links(e->s->parts, e->s->gparts, e->s->sinks, e->s->sparts,
                    e->s->bparts, e->s->nr_parts, e->s->nr_gparts,
                    e->s->nr_sinks, e->s->nr_sparts, e->s->nr_bparts,
                    e->verbose);
#endif

  /* Initial cleaning up session ? */
  if (clean_smoothing_length_values) space_sanitize(e->s);

/* If in parallel, exchange the cell structure, top-level and neighbouring
 * multipoles. To achieve this, free the foreign particle buffers first. */
#ifdef WITH_MPI
  if (e->policy & engine_policy_self_gravity) engine_exchange_top_multipoles(e);

  space_free_foreign_parts(e->s, /*clear_cell_pointers=*/1);

  engine_exchange_cells(e);
#endif

#ifdef SWIFT_DEBUG_CHECKS

  /* Let's check that what we received makes sense */
  if (e->policy & engine_policy_self_gravity) {
    long long counter = 0;

    for (int i = 0; i < e->s->nr_cells; ++i) {
      const struct gravity_tensors *m = &e->s->multipoles_top[i];
      counter += m->m_pole.num_gpart;
    }
    if (counter != e->total_nr_gparts)
      error("Total particles in multipoles inconsistent with engine");
  }
#endif

  /* Re-build the tasks. */
  engine_maketasks(e);

  /* Reallocate freed memory */
#ifdef WITH_MPI
  if (e->free_foreign_when_rebuilding)
    engine_allocate_foreign_particles(e, /*fof=*/0);
#endif

  /* Make the list of top-level cells that have tasks */
  space_list_useful_top_level_cells(e->s);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all cells have been drifted to the current time.
   * That can include cells that have not
   * previously been active on this rank. */
  space_check_drift_point(e->s, e->ti_current,
                          e->policy & engine_policy_self_gravity);

  if (e->policy & engine_policy_self_gravity) {
    for (int k = 0; k < e->s->nr_local_cells; k++)
      cell_check_foreign_multipole(&e->s->cells_top[e->s->local_cells_top[k]]);
  }

  space_check_sort_flags(e->s);

  /* Check whether all the unskip recursion flags are not set */
  space_check_unskip_flags(e->s);
#endif

  /* Run through the tasks and mark as skip or not. */
  if (engine_marktasks(e))
    error("engine_marktasks failed after space_rebuild.");

  /* Print the status of the system */
  if (e->verbose) engine_print_task_counts(e);

  /* Clear the counters of updates since the last rebuild */
  e->updates_since_rebuild = 0;
  e->g_updates_since_rebuild = 0;
  e->s_updates_since_rebuild = 0;
  e->sink_updates_since_rebuild = 0;
  e->b_updates_since_rebuild = 0;

  /* Flag that a rebuild has taken place */
  e->step_props |= engine_step_prop_rebuild;

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Prepare the #engine by re-building the cells and tasks.
 *
 * @param e The #engine to prepare.
 *
 * @return 1 if the function drifted all particles 0 if not.
 */
int engine_prepare(struct engine *e) {

  TIMER_TIC2;
  const ticks tic = getticks();

  int drifted_all = 0;
  int repartitioned = 0;

  /* Unskip active tasks and check for rebuild */
  if (!e->forcerebuild && !e->forcerepart && !e->restarting) engine_unskip(e);

  const ticks tic3 = getticks();

#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &e->forcerebuild, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
#endif

  if (e->verbose)
    message("Communicating rebuild flag took %.3f %s.",
            clocks_from_ticks(getticks() - tic3), clocks_getunit());

  /* Perform FOF search to seed black holes. Only if there is a rebuild coming
   * and no repartitioing. */
  if (e->policy & engine_policy_fof && e->forcerebuild && !e->forcerepart &&
      e->run_fof && e->fof_properties->seed_black_holes_enabled) {

    /* Let's start by drifting everybody to the current time */
    engine_drift_all(e, /*drift_mpole=*/0);
    drifted_all = 1;

    engine_fof(e, e->dump_catalogue_when_seeding, /*dump_debug=*/0,
               /*seed_black_holes=*/1, /*foreign buffers allocated=*/1);

    if (e->dump_catalogue_when_seeding) e->snapshot_output_count++;
  }

  /* Perform particle splitting. Only if there is a rebuild coming and no
     repartitioning. */
  if (!e->restarting && e->forcerebuild && !e->forcerepart && e->step > 1) {

    /* Let's start by drifting everybody to the current time */
    if (!drifted_all) engine_drift_all(e, /*drift_mpole=*/0);
    drifted_all = 1;

    engine_split_gas_particles(e);
  }

  /* Do we need repartitioning ? */
  if (e->forcerepart) {

    /* Let's start by drifting everybody to the current time */
    engine_drift_all(e, /*drift_mpole=*/0);
    drifted_all = 1;

    /* Free the PM grid */
    if ((e->policy & engine_policy_self_gravity) && e->s->periodic)
      pm_mesh_free(e->mesh);

    /* And repartition */
    engine_repartition(e);
    repartitioned = 1;

    /* Reallocate the mesh */
    if ((e->policy & engine_policy_self_gravity) && e->s->periodic)
      pm_mesh_allocate(e->mesh);
  }

  /* Do we need rebuilding ? */
  if (e->forcerebuild) {

    /* Let's start by drifting everybody to the current time */
    if (!e->restarting && !drifted_all) engine_drift_all(e, /*drift_mpole=*/0);

    drifted_all = 1;

    /* And rebuild */
    engine_rebuild(e, repartitioned, 0);
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (e->forcerepart || e->forcerebuild) {
    /* Check that all cells have been drifted to the current time.
     * That can include cells that have not previously been active on this
     * rank. Skip if haven't got any cells (yet). */
    if (e->s->cells_top != NULL)
      space_check_drift_point(e->s, e->ti_current,
                              e->policy & engine_policy_self_gravity);
  }
#endif

  /* Re-rank the tasks every now and then. XXX this never executes. */
  if (e->tasks_age % engine_tasksreweight == 1) {
    scheduler_reweight(&e->sched, e->verbose);
  }
  e->tasks_age += 1;

  TIMER_TOC2(timer_prepare);

  if (e->verbose)
    message("took %.3f %s (including unskip, rebuild and reweight).",
            clocks_from_ticks(getticks() - tic), clocks_getunit());

  return drifted_all;
}

/**
 * @brief Implements a barrier for the #runner threads.
 *
 * @param e The #engine.
 */
void engine_barrier(struct engine *e) {

  /* Wait at the wait barrier. */
  swift_barrier_wait(&e->wait_barrier);

  /* Wait at the run barrier. */
  swift_barrier_wait(&e->run_barrier);
}

/**
 * @brief Print the conserved quantities statistics to a log file
 *
 * @param e The #engine.
 */
void engine_print_stats(struct engine *e) {

  const ticks tic = getticks();

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all cells have been drifted to the current time.
   * That can include cells that have not
   * previously been active on this rank. */
  space_check_drift_point(e->s, e->ti_current, /*chek_mpoles=*/0);

  /* Be verbose about this */
  if (e->nodeID == 0) {
    if (e->policy & engine_policy_cosmology)
      message("Saving statistics at a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Saving statistics at t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }
#else
  if (e->verbose) {
    if (e->policy & engine_policy_cosmology)
      message("Saving statistics at a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Saving statistics at t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }
#endif

  struct statistics stats;
  stats_init(&stats);

  /* Collect the stats on this node */
  stats_collect(e->s, &stats);

/* Aggregate the data from the different nodes. */
#ifdef WITH_MPI
  struct statistics global_stats;
  stats_init(&global_stats);

  if (MPI_Reduce(&stats, &global_stats, 1, statistics_mpi_type,
                 statistics_mpi_reduce_op, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    error("Failed to aggregate stats.");
#else
  struct statistics global_stats = stats;
#endif

  /* Finalize operations */
  stats_finalize(&global_stats);

  /* Print info */
  if (e->nodeID == 0)
    stats_write_to_file(e->file_stats, &global_stats, e->time, e->cosmology->a,
                        e->cosmology->z, e->step);

  /* Flag that we dumped some statistics */
  e->step_props |= engine_step_prop_statistics;

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Sets all the force, drift and kick tasks to be skipped.
 *
 * @param e The #engine to act on.
 */
void engine_skip_force_and_kick(struct engine *e) {

  struct task *tasks = e->sched.tasks;
  const int nr_tasks = e->sched.nr_tasks;

  for (int i = 0; i < nr_tasks; ++i) {

    struct task *t = &tasks[i];

    /* Skip everything that updates the particles */
    if (t->type == task_type_drift_part || t->type == task_type_drift_gpart ||
        t->type == task_type_drift_spart || t->type == task_type_drift_bpart ||
        t->type == task_type_drift_sink || t->type == task_type_kick1 ||
        t->type == task_type_kick2 || t->type == task_type_timestep ||
        t->type == task_type_timestep_limiter ||
        t->type == task_type_timestep_sync || t->type == task_type_collect ||
        t->type == task_type_end_hydro_force || t->type == task_type_cooling ||
        t->type == task_type_stars_in || t->type == task_type_stars_out ||
        t->type == task_type_star_formation ||
        t->type == task_type_star_formation_sink ||
        t->type == task_type_stars_resort || t->type == task_type_extra_ghost ||
        t->type == task_type_stars_ghost ||
        t->type == task_type_stars_ghost_in ||
        t->type == task_type_stars_ghost_out || t->type == task_type_sink_in ||
        t->type == task_type_sink_ghost1 || t->type == task_type_sink_ghost2 ||
        t->type == task_type_sink_formation || t->type == task_type_sink_out ||
        t->type == task_type_stars_prep_ghost1 ||
        t->type == task_type_hydro_prep_ghost1 ||
        t->type == task_type_stars_prep_ghost2 ||
        t->type == task_type_bh_swallow_ghost1 ||
        t->type == task_type_bh_swallow_ghost2 ||
        t->type == task_type_bh_swallow_ghost3 || t->type == task_type_bh_in ||
        t->type == task_type_bh_out || t->type == task_type_rt_ghost1 ||
        t->type == task_type_rt_ghost2 || t->type == task_type_rt_tchem ||
        t->type == task_type_rt_advance_cell_time ||
        t->type == task_type_neutrino_weight || t->type == task_type_csds ||
        t->subtype == task_subtype_force ||
        t->subtype == task_subtype_limiter ||
        t->subtype == task_subtype_gradient ||
        t->subtype == task_subtype_stars_prep1 ||
        t->subtype == task_subtype_stars_prep2 ||
        t->subtype == task_subtype_stars_feedback ||
        t->subtype == task_subtype_bh_feedback ||
        t->subtype == task_subtype_bh_swallow ||
        t->subtype == task_subtype_do_gas_swallow ||
        t->subtype == task_subtype_do_bh_swallow ||
        t->subtype == task_subtype_bpart_rho ||
        t->subtype == task_subtype_part_swallow ||
        t->subtype == task_subtype_bpart_merger ||
        t->subtype == task_subtype_bpart_feedback ||
        t->subtype == task_subtype_sink_swallow ||
        t->subtype == task_subtype_sink_do_sink_swallow ||
        t->subtype == task_subtype_sink_do_gas_swallow ||
        t->subtype == task_subtype_tend || t->subtype == task_subtype_rho ||
        t->subtype == task_subtype_spart_density ||
        t->subtype == task_subtype_part_prep1 ||
        t->subtype == task_subtype_spart_prep2 ||
        t->subtype == task_subtype_sf_counts ||
        t->subtype == task_subtype_rt_gradient ||
        t->subtype == task_subtype_rt_transport)
      t->skip = 1;
  }

  /* Run through the cells and clear some flags. */
  space_map_cells_pre(e->s, 1, cell_clear_drift_flags, NULL);
  space_map_cells_pre(e->s, 1, cell_clear_limiter_flags, NULL);
}

/**
 * @brief Sets all the drift and first kick tasks to be skipped.
 *
 * @param e The #engine to act on.
 */
void engine_skip_drift(struct engine *e) {

  struct task *tasks = e->sched.tasks;
  const int nr_tasks = e->sched.nr_tasks;

  for (int i = 0; i < nr_tasks; ++i) {

    struct task *t = &tasks[i];

    /* Skip everything that moves the particles */
    if (t->type == task_type_drift_part || t->type == task_type_drift_gpart ||
        t->type == task_type_drift_spart || t->type == task_type_drift_bpart ||
        t->type == task_type_drift_sink)
      t->skip = 1;
  }

  /* Run through the cells and clear some flags. */
  space_map_cells_pre(e->s, 1, cell_clear_drift_flags, NULL);
}

/**
 * @brief Launch the runners.
 *
 * @param e The #engine.
 * @param call What kind of tasks are we running? (For time analysis)
 */
void engine_launch(struct engine *e, const char *call) {
  const ticks tic = getticks();

#ifdef SWIFT_DEBUG_CHECKS
  /* Re-set all the cell task counters to 0 */
  space_reset_task_counters(e->s);
#endif

  /* reset the active time counters for the runners */
  for (int i = 0; i < e->nr_threads; ++i) {
    runner_reset_active_time(&e->runners[i]);
  }

  /* Prepare the scheduler. */
  atomic_inc(&e->sched.waiting);

  /* Cry havoc and let loose the dogs of war. */
  swift_barrier_wait(&e->run_barrier);

  /* Load the tasks. */
  scheduler_start(&e->sched);

  /* Remove the safeguard. */
  pthread_mutex_lock(&e->sched.sleep_mutex);
  atomic_dec(&e->sched.waiting);
  pthread_cond_broadcast(&e->sched.sleep_cond);
  pthread_mutex_unlock(&e->sched.sleep_mutex);

  /* Sit back and wait for the runners to come home. */
  swift_barrier_wait(&e->wait_barrier);

  /* Store the wallclock time */
  e->sched.total_ticks += getticks() - tic;

  /* accumulate active counts for all runners */
  ticks active_time = 0;
  for (int i = 0; i < e->nr_threads; ++i) {
    active_time += runner_get_active_time(&e->runners[i]);
  }
  e->sched.deadtime.active_ticks += active_time;
  e->sched.deadtime.waiting_ticks += getticks() - tic;

#ifdef SWIFT_DEBUG_CHECKS
  e->sched.last_successful_task_fetch = 0LL;
#endif

  if (e->verbose)
    message("(%s) took %.3f %s.", call, clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Calls the 'first init' function on the particles of all types.
 *
 * @param e The #engine.
 */
void engine_first_init_particles(struct engine *e) {

  const ticks tic = getticks();

  /* Set the particles in a state where they are ready for a run. */
  space_first_init_parts(e->s, e->verbose);
  space_first_init_gparts(e->s, e->verbose);
  space_first_init_sparts(e->s, e->verbose);
  space_first_init_bparts(e->s, e->verbose);
  space_first_init_sinks(e->s, e->verbose);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Compute the maximal ID of any #part in the run.
 *
 * @param e The #engine.
 */
void engine_get_max_ids(struct engine *e) {

  e->max_parts_id = space_get_max_parts_id(e->s);

#ifdef WITH_MPI
  MPI_Allreduce(MPI_IN_PLACE, &e->max_parts_id, 1, MPI_LONG_LONG_INT, MPI_MAX,
                MPI_COMM_WORLD);
#endif
}

/**
 * @brief Run the radiative transfer sub-cycles outside the
 * regular time-steps.
 *
 * @param e The #engine
 **/
void engine_run_rt_sub_cycles(struct engine *e) {

  /* Do we have work to do? */
  if (!(e->policy & engine_policy_rt)) return;

  /* Note that if running without sub-cycles, no RT-specific timestep data will
   * be written to screen or to the RT subcycles timestep data file. It's
   * meaningless to do so, as all data will already be contained in the normal
   * timesteps file. */
  if (e->max_nr_rt_subcycles <= 1) return;

  /* Get the subcycling step */
  const integertime_t rt_step_size = e->ti_rt_end_min - e->ti_current;
  if (rt_step_size == 0) {
    /* When we arrive at the final step, the rt_step_size can be == 0 */
    if (!engine_is_done(e)) error("Got rt_step_size = 0");
    return;
  }

  /* At this point, the non-RT ti_end_min is up-to-date. Use that and
   * the time of the previous regular step to get how many subcycles
   * we need. */
  const int nr_rt_cycles = (e->ti_end_min - e->ti_current) / rt_step_size;
  /* You can't check here that the number of cycles is exactly the number
   * you fixed it to be. E.g. stars or gravity may reduce the time step
   * sizes for some main steps such that they coincide with the RT bins,
   * yielding effectively no subcycles. (At least for low numbers.) */

  if (nr_rt_cycles < 0) {
    error(
        "Got negative nr of sub-cycles??? ti_rt_end_min = %lld ti_current = "
        "%lld rt_step_size = %lld",
        e->ti_rt_end_min, e->ti_current, rt_step_size);
  } else if (nr_rt_cycles == 0) {
    /* This can happen if in the previous main step no RT/hydro updates
     * happened, but something else (e.g. stars, gravity) only. In this
     * case, exit early. */
    return;
  }

  /* Get some time variables for printouts. Don't update the ones in the
   * engine like in the regular step, or the outputs in the regular steps
   * will be wrong. */
  /* think cosmology one day: needs adapting here */
  /* Also needs adapting further below - we print out current values of a
   * and z. They need to be updated in the engine. */
  if (e->policy & engine_policy_cosmology)
    error("Can't run RT subcycling with cosmology yet");
  const double dt_subcycle = rt_step_size * e->time_base;
  double time = e->ti_current_subcycle * e->time_base + e->time_begin;

  /* Keep track and accumulate the deadtime over all sub-cycles. */
  /* We need to manually put this back in the engine struct when
   * the sub-cycling is completed. */
  double global_deadtime_acc = e->global_deadtime;

  /* Collect and print info before it's gone */
  engine_collect_end_of_sub_cycle(e);

  if (e->nodeID == 0) {

    printf(
        " [rt-sc] %-4d %12e %11.6f %11.6f %13e %4d %4d %12lld %12s %12s "
        "%12s %12s %21s %6s %17s\n",
        0, e->time, e->cosmology->a, e->cosmology->z, dt_subcycle,
        e->min_active_bin_subcycle, e->max_active_bin_subcycle, e->rt_updates,
        /*g, s, sink, bh updates=*/"-", "-", "-", "-", /*wallclock_time=*/"-",
        /*props=*/"-", /*dead_time=*/"-");
#ifdef SWIFT_DEBUG_CHECKS
    fflush(stdout);
#endif

    if (!e->restarting) {
      fprintf(
          e->file_rt_subcycles,
          "  %6d %9d %14e %12.7f %12.7f %14e %4d %4d %12lld %21.3f %17.3f\n",
          e->step, 0, time, e->cosmology->a, e->cosmology->z, dt_subcycle,
          e->min_active_bin_subcycle, e->max_active_bin_subcycle, e->rt_updates,
          /*wall-clock time=*/-1.f, /*deadtime=*/-1.f);
    }
#ifdef SWIFT_DEBUG_CHECKS
    fflush(e->file_rt_subcycles);
#endif
  }

  /* Take note of the (integer) time until which the radiative transfer
   * has been integrated so far. At the start of the sub-cycling, this
   * should be e->ti_current_subcycle + dt_rt_min, since the first (i.e.
   * zeroth) RT cycle has been completed during the regular step.
   * This is used for a consistency/debugging check. */
  integertime_t rt_integration_end = e->ti_current_subcycle + rt_step_size;

  for (int sub_cycle = 1; sub_cycle < nr_rt_cycles; ++sub_cycle) {

    /* Keep track of the wall-clock time of each additional sub-cycle. */
    struct clocks_time time1, time2;
    clocks_gettime(&time1);

    /* reset the deadtime information in the scheduler */
    e->sched.deadtime.active_ticks = 0;
    e->sched.deadtime.waiting_ticks = 0;

    /* Set and re-set times, bins, etc. */
    e->rt_updates = 0ll;
    integertime_t ti_subcycle_old = e->ti_current_subcycle;
    e->ti_current_subcycle = e->ti_current + sub_cycle * rt_step_size;
    e->max_active_bin_subcycle = get_max_active_bin(e->ti_current_subcycle);
    e->min_active_bin_subcycle =
        get_min_active_bin(e->ti_current_subcycle, ti_subcycle_old);
    /* think cosmology one day: needs adapting here */
    if (e->policy & engine_policy_cosmology)
      error("Can't run RT subcycling with cosmology yet");
    time = e->ti_current_subcycle * e->time_base + e->time_begin;

    /* Do the actual work now. */
    engine_unskip_rt_sub_cycle(e);
    TIMER_TIC;
    engine_launch(e, "cycles");
    TIMER_TOC(timer_runners);

    /* Compute the local accumulated deadtime. */
    const ticks deadticks = (e->nr_threads * e->sched.deadtime.waiting_ticks) -
                            e->sched.deadtime.active_ticks;
    e->local_deadtime = clocks_from_ticks(deadticks);

    /* Collect number of updates and print */
    engine_collect_end_of_sub_cycle(e);

    /* Add our sub-cycling deadtime. */
    global_deadtime_acc += e->global_deadtime;

    /* Keep track how far we have integrated over. */
    rt_integration_end += rt_step_size;

    if (e->nodeID == 0) {

      const double dead_time =
          e->global_deadtime / (e->nr_nodes * e->nr_threads);

      /* engine_step() stores the wallclock time in the engine struct.
       * Don't do that here - we want the full step to include the full
       * duration of the step, which includes all sub-cycles. (Also it
       * would be overwritten anyway.) */
      clocks_gettime(&time2);
      const float wallclock_time = (float)clocks_diff(&time1, &time2);

      printf(
          " [rt-sc] %-4d %12e %11.6f %11.6f %13e %4d %4d %12lld %12s %12s "
          "%12s %12s %21.3f %6s %17.3f\n",
          sub_cycle, time, e->cosmology->a, e->cosmology->z, dt_subcycle,
          e->min_active_bin_subcycle, e->max_active_bin_subcycle, e->rt_updates,
          /*g, s, sink, bh updates=*/"-", "-", "-", "-", wallclock_time,
          /*props=*/"-", dead_time);
#ifdef SWIFT_DEBUG_CHECKS
      fflush(stdout);
#endif
      fprintf(
          e->file_rt_subcycles,
          "  %6d %9d %14e %12.7f %12.7f %14e %4d %4d %12lld %21.3f %17.3f\n",
          e->step, sub_cycle, time, e->cosmology->a, e->cosmology->z,
          dt_subcycle, e->min_active_bin_subcycle, e->max_active_bin_subcycle,
          e->rt_updates, wallclock_time, dead_time);
#ifdef SWIFT_DEBUG_CHECKS
      fflush(e->file_rt_subcycles);
#endif
    }
  }

  if (rt_integration_end != e->ti_end_min)
    error(
        "End of sub-cycling doesn't add up: got %lld should have %lld. Started "
        "at ti_current = %lld dt_rt = %lld cycles = %d",
        rt_integration_end, e->ti_end_min, e->ti_current, rt_step_size,
        nr_rt_cycles);

  /* Once we're done, clean up after ourselves */
  e->rt_updates = 0ll;
  e->global_deadtime = global_deadtime_acc;
}

/**
 * @brief Initialises the particles and set them in a state ready to move
 *forward in time.
 *
 * @param e The #engine
 * @param flag_entropy_ICs Did the 'Internal Energy' of the particles actually
 * contain entropy ?
 * @param clean_h_values Are we cleaning up the values of h before building
 */
void engine_init_particles(struct engine *e, int flag_entropy_ICs,
                           int clean_h_values) {

  struct space *s = e->s;

  struct clocks_time time1, time2;
  clocks_gettime(&time1);

  /* reset the deadtime information in the scheduler */
  e->sched.deadtime.active_ticks = 0;
  e->sched.deadtime.waiting_ticks = 0;

  /* Update the softening lengths */
  if (e->policy & engine_policy_self_gravity)
    gravity_props_update(e->gravity_properties, e->cosmology);

  /* Udpate the hydro properties */
  if (e->policy & engine_policy_hydro)
    hydro_props_update(e->hydro_properties, e->gravity_properties,
                       e->cosmology);

  /* Start by setting the particles in a good state */
  if (e->nodeID == 0) message("Setting particles to a valid state...");
  engine_first_init_particles(e);

  if (e->nodeID == 0)
    message("Computing initial gas densities and approximate gravity.");

  /* Construct all cells and tasks to start everything */
  engine_rebuild(e, 0, clean_h_values);

  /* Compute the mesh forces for the first time */
  if ((e->policy & engine_policy_self_gravity) && e->s->periodic) {

    /* Compute mesh forces */
    pm_mesh_compute_potential(e->mesh, e->s, &e->threadpool, e->verbose);

    /* Compute mesh time-step length */
    engine_recompute_displacement_constraint(e);

    e->step_props |= engine_step_prop_mesh;
  }

  /* No time integration. We just want the density and ghosts */
  engine_skip_force_and_kick(e);

  /* Print the number of active tasks ? */
  if (e->verbose) engine_print_task_counts(e);

  /* Init the particle data (by hand). */
  space_init_parts(s, e->verbose);
  space_init_gparts(s, e->verbose);
  space_init_sparts(s, e->verbose);
  space_init_bparts(s, e->verbose);
  space_init_sinks(s, e->verbose);

  /* Update the cooling function */
  if ((e->policy & engine_policy_cooling) ||
      (e->policy & engine_policy_temperature))
    cooling_update(e->physical_constants, e->cosmology, e->pressure_floor_props,
                   e->cooling_func, e->s, e->time);

#ifdef WITH_CSDS
  if (e->policy & engine_policy_csds) {
    /* Mark the first time step in the particle csds file. */
    if (e->policy & engine_policy_cosmology) {
      csds_log_timestamp(e->csds, e->ti_current, e->cosmology->a,
                         &e->csds->timestamp_offset);
    } else {
      csds_log_timestamp(e->csds, e->ti_current, e->time,
                         &e->csds->timestamp_offset);
    }
    /* Make sure that we have enough space in the particle csds file
     * to store the particles in current time step. */
    csds_ensure_size(e->csds, e);
    csds_write_description(e->csds, e);
  }
#endif

  /* Now, launch the calculation */
  TIMER_TIC;
  engine_launch(e, "tasks");
  TIMER_TOC(timer_runners);

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /* Run the brute-force hydro calculation for some parts */
  if (e->policy & engine_policy_hydro)
    hydro_exact_density_compute(e->s, e, /*check_force=*/0);

  /* Check the accuracy of the hydro calculation */
  if (e->policy & engine_policy_hydro)
    hydro_exact_density_check(e->s, e, /*rel_tol=*/1e-3, /*check_force=*/0);
#endif

  /* Apply some conversions (e.g. internal energy -> entropy) */
  if (!flag_entropy_ICs) {

    if (e->nodeID == 0) message("Converting internal energy variable.");

    space_convert_quantities(e->s, e->verbose);

    /* Correct what we did (e.g. in PE-SPH, need to recompute rho_bar) */
    if (hydro_need_extra_init_loop) {
      engine_marktasks(e);
      engine_skip_force_and_kick(e);
      engine_launch(e, "tasks");
    }
  }

  /* Do some post initialisations */
  space_post_init_parts(e->s, e->verbose);

  /* Apply some RT conversions (e.g. energy -> energy density) */
  if (e->policy & engine_policy_rt)
    space_convert_rt_quantities(e->s, e->verbose);

  /* Collect initial mean mass of each particle type */
  space_collect_mean_masses(e->s, e->verbose);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we have the correct total mass in the top-level multipoles */
  long long num_gpart_mpole = 0;
  if (e->policy & engine_policy_self_gravity) {
    for (int i = 0; i < e->s->nr_cells; ++i)
      num_gpart_mpole += e->s->cells_top[i].grav.multipole->m_pole.num_gpart;
    if (num_gpart_mpole != e->total_nr_gparts)
      error(
          "Top-level multipoles don't contain the total number of gpart "
          "s->nr_gpart=%lld, "
          "m_poles=%lld",
          e->total_nr_gparts, num_gpart_mpole);
  }
#endif

  /* Now time to get ready for the first time-step */
  if (e->nodeID == 0) message("Running initial fake time-step.");

  /* Update the MAC strategy if necessary */
  if (e->policy & engine_policy_self_gravity)
    gravity_props_update_MAC_choice(e->gravity_properties);

  /* Construct all cells again for a new round (need to update h_max) */
  engine_rebuild(e, 0, 0);

  /* No drift this time */
  engine_skip_drift(e);

  /* Init the particle data (by hand). */
  space_init_parts(e->s, e->verbose);
  space_init_gparts(e->s, e->verbose);
  space_init_sparts(e->s, e->verbose);
  space_init_bparts(e->s, e->verbose);
  space_init_sinks(e->s, e->verbose);

  /* Print the number of active tasks ? */
  if (e->verbose) engine_print_task_counts(e);

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Run the brute-force gravity calculation for some gparts */
  if (e->policy & engine_policy_self_gravity)
    gravity_exact_force_compute(e->s, e);
#endif

  scheduler_write_dependencies(&e->sched, e->verbose, e->step);
  scheduler_write_cell_dependencies(&e->sched, e->verbose, e->step);
  if (e->nodeID == 0) scheduler_write_task_level(&e->sched, e->step);

  /* Run the 0th time-step */
  TIMER_TIC2;
  engine_launch(e, "tasks");
  TIMER_TOC2(timer_runners);

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /* Run the brute-force hydro calculation for some parts */
  if (e->policy & engine_policy_hydro)
    hydro_exact_density_compute(e->s, e, /*check_force=*/1);

  /* Check the accuracy of the hydro calculation */
  if (e->policy & engine_policy_hydro)
    hydro_exact_density_check(e->s, e, /*rel_tol=*/1e-3, /*check_force=*/1);
#endif

#ifdef SWIFT_STARS_DENSITY_CHECKS
  /* Run the brute-force stars calculation for some parts */
  if (e->policy & engine_policy_stars) stars_exact_density_compute(e->s, e);

  /* Check the accuracy of the stars calculation */
  if (e->policy & engine_policy_stars)
    stars_exact_density_check(e->s, e, /*rel_tol=*/1e-3);
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Check the accuracy of the gravity calculation */
  if (e->policy & engine_policy_self_gravity)
    gravity_exact_force_check(e->s, e, 1e-1);
#endif

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure all woken-up particles have been processed */
  space_check_limiter(e->s);
  space_check_swallow(e->s);
#endif

  /* Compute the local accumulated deadtime. */
  const ticks deadticks = (e->nr_threads * e->sched.deadtime.waiting_ticks) -
                          e->sched.deadtime.active_ticks;
  e->local_deadtime = clocks_from_ticks(deadticks);

  /* Recover the (integer) end of the next time-step */
  engine_collect_end_of_step(e, 1);

  /* Check if any particles have the same position. This is not
   * allowed (/0) so we abort.*/
  if (s->nr_parts > 0) {

    /* Sorting should put the same positions next to each other... */
    int failed = 0;
    double *prev_x = s->parts[0].x;
    long long *prev_id = &s->parts[0].id;
    for (size_t k = 1; k < s->nr_parts; k++) {

      /* Ignore fake buffer particles for on-the-fly creation */
      if (s->parts[k].time_bin == time_bin_not_created) continue;

      if (prev_x[0] == s->parts[k].x[0] && prev_x[1] == s->parts[k].x[1] &&
          prev_x[2] == s->parts[k].x[2]) {
        if (e->verbose)
          message("Two particles occupy location: %f %f %f id=%lld id=%lld",
                  prev_x[0], prev_x[1], prev_x[2], *prev_id, s->parts[k].id);
        failed++;
      }
      prev_x = s->parts[k].x;
      prev_id = &s->parts[k].id;
    }
    if (failed > 0)
      error(
          "Have %d particle pairs with the same locations.\n"
          "Cannot continue",
          failed);
  }

  /* Also check any gparts. This is not supposed to be fatal so only warn. */
  if (s->nr_gparts > 0) {
    int failed = 0;
    double *prev_x = s->gparts[0].x;
    for (size_t k = 1; k < s->nr_gparts; k++) {

      /* Ignore fake buffer particles for on-the-fly creation */
      if (s->gparts[k].time_bin == time_bin_not_created) continue;

      if (prev_x[0] == s->gparts[k].x[0] && prev_x[1] == s->gparts[k].x[1] &&
          prev_x[2] == s->gparts[k].x[2]) {
        if (e->verbose)
          message("Two gparts occupy location: %f %f %f / %f %f %f", prev_x[0],
                  prev_x[1], prev_x[2], s->gparts[k].x[0], s->gparts[k].x[1],
                  s->gparts[k].x[2]);
        failed++;
      }
      prev_x = s->gparts[k].x;
    }
    if (failed > 0)
      message(
          "WARNING: found %d gpart pairs at the same location. "
          "That is not optimal",
          failed);
  }

  /* Check the top-level cell h_max matches the particles as these can be
   * updated in the the ghost tasks (only a problem if the ICs estimates for h
   * are too small). Note this must be followed by a rebuild as sub-cells will
   * not be updated until that is done. */
  if (s->cells_top != NULL && s->nr_parts > 0) {
    for (int i = 0; i < s->nr_cells; i++) {
      struct cell *c = &s->cells_top[i];
      if (c->nodeID == engine_rank && c->hydro.count > 0) {
        float part_h_max = c->hydro.parts[0].h;
        for (int k = 1; k < c->hydro.count; k++) {
          if (c->hydro.parts[k].h > part_h_max)
            part_h_max = c->hydro.parts[k].h;
        }
        c->hydro.h_max = max(part_h_max, c->hydro.h_max);
      }
    }
  }

  if (s->cells_top != NULL && s->nr_sparts > 0) {
    for (int i = 0; i < s->nr_cells; i++) {
      struct cell *c = &s->cells_top[i];
      if (c->nodeID == engine_rank && c->stars.count > 0) {
        float spart_h_max = c->stars.parts[0].h;
        for (int k = 1; k < c->stars.count; k++) {
          if (c->stars.parts[k].h > spart_h_max)
            spart_h_max = c->stars.parts[k].h;
        }
        c->stars.h_max = max(spart_h_max, c->stars.h_max);
      }
    }
  }

  if (s->cells_top != NULL && s->nr_sinks > 0) {
    for (int i = 0; i < s->nr_cells; i++) {
      struct cell *c = &s->cells_top[i];
      if (c->nodeID == engine_rank && c->sinks.count > 0) {
        float sink_h_max = c->sinks.parts[0].r_cut;
        for (int k = 1; k < c->sinks.count; k++) {
          if (c->sinks.parts[k].r_cut > sink_h_max)
            sink_h_max = c->sinks.parts[k].r_cut;
        }
        c->sinks.r_cut_max = max(sink_h_max, c->sinks.r_cut_max);
      }
    }
  }

  /* Run the RT sub-cycles now. */
  engine_run_rt_sub_cycles(e);

  clocks_gettime(&time2);

#ifdef SWIFT_DEBUG_CHECKS
  space_check_timesteps(e->s);
  part_verify_links(e->s->parts, e->s->gparts, e->s->sinks, e->s->sparts,
                    e->s->bparts, e->s->nr_parts, e->s->nr_gparts,
                    e->s->nr_sinks, e->s->nr_sparts, e->s->nr_bparts,
                    e->verbose);
#endif

  /* Gather the max IDs at this stage */
  engine_get_max_ids(e);

  /* Ready to go */
  e->step = 0;
  e->forcerebuild = 1;
  e->wallclock_time = (float)clocks_diff(&time1, &time2);

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  e->force_checks_snapshot_flag = 0;
#endif

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* if we're running the debug RT scheme, do some checks after every step,
   * and reset debugging flags now. */
  if (e->policy & engine_policy_rt) {
    rt_debugging_checks_end_of_step(e, e->verbose);
  }
#endif

  if (e->verbose) message("took %.3f %s.", e->wallclock_time, clocks_getunit());
}

/**
 * @brief Let the #engine loose to compute the forces.
 *
 * @param e The #engine.
 * @return Should the run stop after this step?
 */
int engine_step(struct engine *e) {

  TIMER_TIC2;

  struct clocks_time time1, time2;
  clocks_gettime(&time1);

  /* reset the deadtime information in the scheduler */
  e->sched.deadtime.active_ticks = 0;
  e->sched.deadtime.waiting_ticks = 0;

#if defined(SWIFT_MPIUSE_REPORTS) && defined(WITH_MPI)
  /* We may want to compare times across ranks, so make sure all steps start
   * at the same time, just different ticks. */
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  e->tic_step = getticks();

  if (e->nodeID == 0) {

    const double dead_time = e->global_deadtime / (e->nr_nodes * e->nr_threads);

    const ticks tic_files = getticks();

    /* Print some information to the screen */
    printf(
        "  %6d %14e %12.7f %12.7f %14e %4d %4d %12lld %12lld %12lld "
        "%12lld %12lld %21.3f %6d %17.3f\n",
        e->step, e->time, e->cosmology->a, e->cosmology->z, e->time_step,
        e->min_active_bin, e->max_active_bin, e->updates, e->g_updates,
        e->s_updates, e->sink_updates, e->b_updates, e->wallclock_time,
        e->step_props, dead_time);
#ifdef SWIFT_DEBUG_CHECKS
    fflush(stdout);
#endif

    /* Write the star formation information to the file */
    if (e->policy & engine_policy_star_formation) {

      star_formation_logger_write_to_log_file(e->sfh_logger, e->time,
                                              e->cosmology->a, e->cosmology->z,
                                              e->sfh, e->step);

#ifdef SWIFT_DEBUG_CHECKS
      fflush(e->sfh_logger);
#else
      if (e->step % 32 == 0) fflush(e->sfh_logger);
#endif
    }

    if (!e->restarting)
      fprintf(
          e->file_timesteps,
          "  %6d %14e %12.7f %12.7f %14e %4d %4d %12lld %12lld %12lld %12lld "
          "%12lld %21.3f %6d %17.3f\n",
          e->step, e->time, e->cosmology->a, e->cosmology->z, e->time_step,
          e->min_active_bin, e->max_active_bin, e->updates, e->g_updates,
          e->s_updates, e->sink_updates, e->b_updates, e->wallclock_time,
          e->step_props, dead_time);
#ifdef SWIFT_DEBUG_CHECKS
    fflush(e->file_timesteps);
#endif

    if (e->verbose)
      message("Writing step info to files took %.3f %s",
              clocks_from_ticks(getticks() - tic_files), clocks_getunit());
  }

  /* When restarting, we may have had some i/o to do on the step
   * where we decided to stop. We have to do this now.
   * We need some cells to exist but not the whole task stuff. */
  if (e->restarting) space_rebuild(e->s, 0, e->verbose);
  if (e->restarting) engine_io(e);

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* If we're restarting, clean up some flags and counters first. If would
   * usually be done at the end of the step, but the restart dump
   * interrupts it. */
  if (e->restarting && (e->policy & engine_policy_rt))
    rt_debugging_checks_end_of_step(e, e->verbose);
#endif

  /* Move forward in time */
  e->ti_old = e->ti_current;
  e->ti_current = e->ti_end_min;
  e->max_active_bin = get_max_active_bin(e->ti_end_min);
  e->min_active_bin = get_min_active_bin(e->ti_current, e->ti_old);
  e->step += 1;
  engine_current_step = e->step;
  e->step_props = engine_step_prop_none;

  /* RT sub-cycling related time updates */
  e->max_active_bin_subcycle = get_max_active_bin(e->ti_end_min);
  e->min_active_bin_subcycle =
      get_min_active_bin(e->ti_end_min, e->ti_current_subcycle);
  e->ti_current_subcycle = e->ti_end_min;

  /* When restarting, move everyone to the current time. */
  if (e->restarting) engine_drift_all(e, /*drift_mpole=*/1);

  /* Get the physical value of the time and time-step size */
  if (e->policy & engine_policy_cosmology) {
    e->time_old = e->time;
    cosmology_update(e->cosmology, e->physical_constants, e->ti_current);
    e->time = e->cosmology->time;
    e->time_step = e->time - e->time_old;
  } else {
    e->time = e->ti_current * e->time_base + e->time_begin;
    e->time_old = e->ti_old * e->time_base + e->time_begin;
    e->time_step = (e->ti_current - e->ti_old) * e->time_base;
  }

#ifdef WITH_LIGHTCONE
  /* Determine which periodic replications could contribute to the lightcone
     during this time step */
  lightcone_array_prepare_for_step(e->lightcone_array_properties, e->cosmology,
                                   e->ti_earliest_undrifted, e->ti_current);
#endif

  /*****************************************************/
  /* OK, we now know what the next end of time-step is */
  /*****************************************************/

  const ticks tic_updates = getticks();

  /* Update the cooling function */
  if ((e->policy & engine_policy_cooling) ||
      (e->policy & engine_policy_temperature))
    cooling_update(e->physical_constants, e->cosmology, e->pressure_floor_props,
                   e->cooling_func, e->s, e->time);

  /* Update the softening lengths */
  if (e->policy & engine_policy_self_gravity)
    gravity_props_update(e->gravity_properties, e->cosmology);

  /* Udpate the hydro properties */
  if (e->policy & engine_policy_hydro)
    hydro_props_update(e->hydro_properties, e->gravity_properties,
                       e->cosmology);

  /* Check for any snapshot triggers */
  engine_io_check_snapshot_triggers(e);

  if (e->verbose)
    message("Updating general quantities took %.3f %s",
            clocks_from_ticks(getticks() - tic_updates), clocks_getunit());

  /* Trigger a tree-rebuild if we passed the frequency threshold */
  if ((e->policy & engine_policy_self_gravity) &&
      ((double)e->g_updates_since_rebuild >
       ((double)e->total_nr_gparts) * e->gravity_properties->rebuild_frequency))
    e->forcerebuild = 1;

  /* Trigger a FOF black hole seeding? */
  if (e->policy & engine_policy_fof && !e->restarting) {
    if (e->ti_end_min > e->ti_next_fof && e->ti_next_fof > 0) {
      e->run_fof = 1;
      e->forcerebuild = 1;
    }
  }

  /* Trigger a rebuild if we reached a gravity mesh step? */
  if ((e->policy & engine_policy_self_gravity) && e->s->periodic &&
      e->mesh->ti_end_mesh_next == e->ti_current)
    e->forcerebuild = 1;

  /* Do we want a snapshot that will trigger a FOF call? */
  if (e->ti_current + (e->ti_current - e->ti_old) > e->ti_next_snapshot &&
      e->ti_next_snapshot > 0) {
    if ((e->policy & engine_policy_fof) && e->snapshot_invoke_fof) {
      e->forcerebuild = 1;
    }
  }

  /* Trigger a tree-rebuild if the fraction of active gparts is large enough */
  if ((e->policy & engine_policy_self_gravity) && !e->forcerebuild &&
      e->gravity_properties->rebuild_active_fraction <= 1.0f) {

    ticks tic = getticks();

    /* Count the number of active particles */
    size_t nr_gparts = e->s->nr_gparts;
    size_t nr_active_gparts = 0;
    for (size_t i = 0; i < nr_gparts; ++i) {
      struct gpart *gp = &e->s->gparts[i];
      if (gpart_is_active(gp, e)) nr_active_gparts++;
    }

    long long total_nr_active_gparts = nr_active_gparts;
#ifdef WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, &total_nr_active_gparts, 1, MPI_LONG_LONG_INT,
                  MPI_SUM, MPI_COMM_WORLD);
#endif

    if (e->verbose)
      message("Counting active gparts took %.3f %s.",
              clocks_from_ticks(getticks() - tic), clocks_getunit());

    /* Trigger the tree-rebuild? */
    if (((double)total_nr_active_gparts >
         ((double)e->total_nr_gparts) *
             e->gravity_properties->rebuild_active_fraction))
      e->forcerebuild = 1;
  }

#ifdef WITH_CSDS
  if (e->policy & engine_policy_csds) {
    /* Mark the current time step in the particle csds file. */
    if (e->policy & engine_policy_cosmology) {
      csds_log_timestamp(e->csds, e->ti_current, e->cosmology->a,
                         &e->csds->timestamp_offset);
    } else {
      csds_log_timestamp(e->csds, e->ti_current, e->time,
                         &e->csds->timestamp_offset);
    }
    /* Make sure that we have enough space in the particle csds file
     * to store the particles in current time step. */
    csds_ensure_size(e->csds, e);
  }
#endif

  /* Are we drifting everything (a la Gadget/GIZMO) ? */
  if (e->policy & engine_policy_drift_all && !e->forcerebuild)
    engine_drift_all(e, /*drift_mpole=*/1);

  /* Are we reconstructing the multipoles or drifting them ?*/
  if ((e->policy & engine_policy_self_gravity) && !e->forcerebuild) {

    if (e->policy & engine_policy_reconstruct_mpoles)
      engine_reconstruct_multipoles(e);
    else
      engine_drift_top_multipoles(e);
  }

#ifdef WITH_MPI
  /* Repartition the space amongst the nodes? */
  engine_repartition_trigger(e);
#endif

  /* Prepare the tasks to be launched, rebuild or repartition if needed. */
  const int drifted_all = engine_prepare(e);

  /* Dump local cells and active particle counts. */
  // dumpCells("cells", 1, 0, 0, 0, e->s, e->nodeID, e->step);

#ifdef SWIFT_DEBUG_CHECKS
  /* Print the number of active tasks */
  if (e->verbose) engine_print_task_counts(e);

  /* Check that we have the correct total mass in the top-level multipoles */
  long long num_gpart_mpole = 0;
  if (e->policy & engine_policy_self_gravity) {
    for (int i = 0; i < e->s->nr_cells; ++i)
      num_gpart_mpole += e->s->cells_top[i].grav.multipole->m_pole.num_gpart;
    if (num_gpart_mpole != e->total_nr_gparts)
      error(
          "Multipoles don't contain the total number of gpart mpoles=%lld "
          "ngparts=%lld",
          num_gpart_mpole, e->total_nr_gparts);
  }
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Do we need to check if all gparts are active? */
  if (e->force_checks_only_all_active) {
    size_t nr_gparts = e->s->nr_gparts;
    e->all_gparts_active = 1;

    /* Look for inactive gparts */
    for (size_t i = 0; i < nr_gparts; ++i) {
      struct gpart *gp = &e->s->gparts[i];

      /* If one gpart is inactive we can stop. */
      if (!gpart_is_active(gp, e)) {
        e->all_gparts_active = 0;
        break;
      }
    }
  }

  /* Check if we want to run force checks this timestep. */
  if (e->policy & engine_policy_self_gravity) {
    /* Are all gparts active (and the option is selected)? */
    if ((e->all_gparts_active && e->force_checks_only_all_active) ||
        !e->force_checks_only_all_active) {
      /* Is this a snapshot timestep (and the option is selected)? */
      if ((e->force_checks_snapshot_flag &&
           e->force_checks_only_at_snapshots) ||
          !e->force_checks_only_at_snapshots) {

        /* Do checks */
        gravity_exact_force_compute(e->s, e);
      }
    }
  }
#endif

  /* Re-compute the mesh forces? */
  if ((e->policy & engine_policy_self_gravity) && e->s->periodic &&
      e->mesh->ti_end_mesh_next == e->ti_current) {

    /* We might need to drift things */
    if (!drifted_all) engine_drift_all(e, /*drift_mpole=*/0);

    /* ... and recompute */
    pm_mesh_compute_potential(e->mesh, e->s, &e->threadpool, e->verbose);

    /* Check whether we need to update the mesh time-step length */
    engine_recompute_displacement_constraint(e);

    e->step_props |= engine_step_prop_mesh;
  }

  /* Get current CPU times.*/
#ifdef WITH_MPI
  double start_usertime = 0.0;
  double start_systime = 0.0;
  clocks_get_cputimes_used(&start_usertime, &start_systime);
#endif

  /* Write the dependencies */
  if (e->sched.frequency_dependency != 0 &&
      e->step % e->sched.frequency_dependency == 0) {
    scheduler_write_dependencies(&e->sched, e->verbose, e->step);
    scheduler_write_cell_dependencies(&e->sched, e->verbose, e->step);
  }

  /* Write the task levels */
  if (e->sched.frequency_task_levels != 0 &&
      e->step % e->sched.frequency_task_levels == 0)
    scheduler_write_task_level(&e->sched, e->step);

  /* we have to reset the ghost histograms here and not in engine_launch,
     because engine_launch is re-used for the limiter and sync (and we don't
     want to lose the data from the tasks) */
  space_reset_ghost_histograms(e->s);

  /* Start all the tasks. */
  TIMER_TIC;
  engine_launch(e, "tasks");
  TIMER_TOC(timer_runners);

  /* Now record the CPU times used by the tasks. */
#ifdef WITH_MPI
  double end_usertime = 0.0;
  double end_systime = 0.0;
  clocks_get_cputimes_used(&end_usertime, &end_systime);
  e->usertime_last_step = end_usertime - start_usertime;
  e->systime_last_step = end_systime - start_systime;
#endif

#ifdef SWIFT_HYDRO_DENSITY_CHECKS
  /* Run the brute-force hydro calculation for some parts */
  if (e->policy & engine_policy_hydro)
    hydro_exact_density_compute(e->s, e, /*check_force=*/1);

  /* Check the accuracy of the hydro calculation */
  if (e->policy & engine_policy_hydro)
    hydro_exact_density_check(e->s, e, /*rel_tol=*/1e-3, /*check_force=*/1);
#endif

#ifdef SWIFT_STARS_DENSITY_CHECKS
  /* Run the brute-force stars calculation for some parts */
  if (e->policy & engine_policy_stars) stars_exact_density_compute(e->s, e);

  /* Check the accuracy of the stars calculation */
  if (e->policy & engine_policy_stars)
    stars_exact_density_check(e->s, e, /*rel_tol=*/1e-2);
#endif

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  /* Check if we want to run force checks this timestep. */
  if (e->policy & engine_policy_self_gravity) {
    /* Are all gparts active (and the option is selected)? */
    if ((e->all_gparts_active && e->force_checks_only_all_active) ||
        !e->force_checks_only_all_active) {
      /* Is this a snapshot timestep (and the option is selected)? */
      if ((e->force_checks_snapshot_flag &&
           e->force_checks_only_at_snapshots) ||
          !e->force_checks_only_at_snapshots) {

        /* Do checks */
        gravity_exact_force_check(e->s, e, 1e-1);

        /* Reset flag waiting for next output time */
        e->force_checks_snapshot_flag = 0;
      }
    }
  }
#endif

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure all woken-up particles have been processed */
  space_check_limiter(e->s);
  space_check_sort_flags(e->s);
  space_check_swallow(e->s);

  /* Verify that all the unskip flags for the gravity have been cleaned */
  space_check_unskip_flags(e->s);
#endif

  /* Compute the local accumulated deadtime. */
  const ticks deadticks = (e->nr_threads * e->sched.deadtime.waiting_ticks) -
                          e->sched.deadtime.active_ticks;
  e->local_deadtime = clocks_from_ticks(deadticks);

  /* Collect information about the next time-step */
  engine_collect_end_of_step(e, 1);
  e->forcerebuild = e->collect_group1.forcerebuild;
  e->updates_since_rebuild += e->collect_group1.updated;
  e->g_updates_since_rebuild += e->collect_group1.g_updated;
  e->s_updates_since_rebuild += e->collect_group1.s_updated;
  e->sink_updates_since_rebuild += e->collect_group1.sink_updated;
  e->b_updates_since_rebuild += e->collect_group1.b_updated;

  /* Check if we updated all of the particles on this step */
  if ((e->collect_group1.updated == e->total_nr_parts) &&
      (e->collect_group1.g_updated == e->total_nr_gparts) &&
      (e->collect_group1.s_updated == e->total_nr_sparts) &&
      (e->collect_group1.sink_updated == e->total_nr_sinks) &&
      (e->collect_group1.b_updated == e->total_nr_bparts))
    e->ti_earliest_undrifted = e->ti_current;

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify that all cells have correct time-step information */
  space_check_timesteps(e->s);

  if (e->ti_end_min == e->ti_current && e->ti_end_min < max_nr_timesteps)
    error("Obtained a time-step of size 0");
#endif

  /* Run the RT sub-cycling now. */
  engine_run_rt_sub_cycles(e);

#ifdef WITH_CSDS
  if (e->policy & engine_policy_csds && e->verbose)
    message("The CSDS currently uses %f GB of storage",
            e->collect_group1.csds_file_size_gb);
#endif

    /********************************************************/
    /* OK, we are done with the regular stuff. Time for i/o */
    /********************************************************/

#ifdef WITH_LIGHTCONE
  /* Flush lightcone buffers if necessary */
  const int flush = e->flush_lightcone_maps;
  lightcone_array_flush(e->lightcone_array_properties, &(e->threadpool),
                        e->cosmology, e->internal_units, e->snapshot_units,
                        /*flush_map_updates=*/flush, /*flush_particles=*/0,
                        /*end_file=*/0, /*dump_all_shells=*/0);
#endif

  /* Create a restart file if needed. */
  const int force_stop =
      engine_dump_restarts(e, 0, e->restart_onexit && engine_is_done(e));

  /* Is there any form of i/o this step?
   *
   * Note that if the run was forced to stop, we do not dump,
   * we will do so when the run is restarted*/
  if (!force_stop) engine_io(e);

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* if we're running the debug RT scheme, do some checks after every step.
   * Do this after the output so we can safely reset debugging checks now. */
  if (e->policy & engine_policy_rt)
    rt_debugging_checks_end_of_step(e, e->verbose);
#endif

  TIMER_TOC2(timer_step);

  clocks_gettime(&time2);
  e->wallclock_time = (float)clocks_diff(&time1, &time2);

  /* Time in ticks at the end of this step. */
  e->toc_step = getticks();

  return force_stop;
}

/**
 * @brief Returns 1 if the simulation has reached its end point, 0 otherwise
 */
int engine_is_done(struct engine *e) {
  return !(e->ti_current < max_nr_timesteps);
}

void engine_do_reconstruct_multipoles_mapper(void *map_data, int num_elements,
                                             void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct cell *cells = (struct cell *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &cells[ind];
    if (c != NULL && c->nodeID == e->nodeID) {

      /* Construct the multipoles in this cell hierarchy */
      cell_make_multipoles(c, e->ti_current, e->gravity_properties);
    }
  }
}

/**
 * @brief Reconstruct all the multipoles at all the levels in the tree.
 *
 * @param e The #engine.
 */
void engine_reconstruct_multipoles(struct engine *e) {

  const ticks tic = getticks();

#ifdef SWIFT_DEBUG_CHECKS
  if (e->nodeID == 0) {
    if (e->policy & engine_policy_cosmology)
      message("Reconstructing multipoles at a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Reconstructing multipoles at t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }
#endif

  threadpool_map(&e->threadpool, engine_do_reconstruct_multipoles_mapper,
                 e->s->cells_top, e->s->nr_cells, sizeof(struct cell),
                 threadpool_auto_chunk_size, e);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Split the underlying space into regions, construct proxies and
 * distribute the particles where they belong.
 *
 * @param e The #engine.
 * @param initial_partition structure defining the cell partition technique
 */
void engine_split(struct engine *e, struct partition *initial_partition) {

#ifdef WITH_MPI
  const ticks tic = getticks();

  struct space *s = e->s;

  /* Do the initial partition of the cells. */
  partition_initial_partition(initial_partition, e->nodeID, e->nr_nodes, s);

  /* Make the proxies. */
  engine_makeproxies(e);

  /* Turn off the csds to avoid writing the communications to
   * the CSDS (since we haven't properly started the run yet) */
  const int with_csds = e->policy & engine_policy_csds;
  if (with_csds) e->policy &= ~engine_policy_csds;

  /* Move the particles to the ranks they belong to */
  engine_redistribute(e);

  /* Turn it back on */
  if (with_csds) e->policy |= engine_policy_csds;

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

#ifdef DEBUG_INTERACTIONS_STARS
/**
 * @brief Exchange the feedback counters between stars
 * @param e The #engine.
 */
void engine_collect_stars_counter(struct engine *e) {

#ifdef WITH_MPI
  if (e->total_nr_sparts > 1e5) {
    message("WARNING: too many sparts, skipping exchange of counters");
    return;
  }

  /* Get number of sparticles for each rank */
  size_t *n_sparts = (size_t *)malloc(e->nr_nodes * sizeof(size_t));

  int err = MPI_Allgather(&e->s->nr_sparts_foreign, 1, MPI_UNSIGNED_LONG,
                          n_sparts, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) error("Communication failed");

  /* Compute derivated quantities */
  int total = 0;
  int *n_sparts_int = (int *)malloc(e->nr_nodes * sizeof(int));
  int *displs = (int *)malloc(e->nr_nodes * sizeof(int));
  for (int i = 0; i < e->nr_nodes; i++) {
    displs[i] = total;
    total += n_sparts[i];
    n_sparts_int[i] = n_sparts[i];
  }

  /* Get all sparticles */
  struct spart *sparts =
      (struct spart *)swift_malloc("sparts", total * sizeof(struct spart));
  err = MPI_Allgatherv(e->s->sparts_foreign, e->s->nr_sparts_foreign,
                       spart_mpi_type, sparts, n_sparts_int, displs,
                       spart_mpi_type, MPI_COMM_WORLD);
  if (err != MPI_SUCCESS) error("Communication failed");

  /* Reset counters */
  for (size_t i = 0; i < e->s->nr_sparts_foreign; i++) {
    e->s->sparts_foreign[i].num_ngb_feedback = 0;
  }

  /* Update counters */
  struct spart *local_sparts = e->s->sparts;
  for (size_t i = 0; i < e->s->nr_sparts; i++) {
    const long long id_i = local_sparts[i].id;

    for (int j = 0; j < total; j++) {
      const long long id_j = sparts[j].id;

      if (id_j == id_i) {
        if (j >= displs[engine_rank] &&
            j < displs[engine_rank] + n_sparts_int[engine_rank]) {
          error(
              "Found a local spart in foreign cell ID=%lli: j=%i, displs=%i, "
              "n_sparts=%i",
              id_j, j, displs[engine_rank], n_sparts_int[engine_rank]);
        }

        local_sparts[i].num_ngb_feedback += sparts[j].num_ngb_feedback;
      }
    }
  }

  free(n_sparts);
  free(n_sparts_int);
  swift_free("sparts", sparts);
#endif
}

#endif

#ifdef HAVE_SETAFFINITY
/**
 * @brief Returns the initial affinity the main thread is using.
 */
cpu_set_t *engine_entry_affinity(void) {

  static int use_entry_affinity = 0;
  static cpu_set_t entry_affinity;

  if (!use_entry_affinity) {
    pthread_t engine = pthread_self();
    pthread_getaffinity_np(engine, sizeof(entry_affinity), &entry_affinity);
    use_entry_affinity = 1;
  }

  return &entry_affinity;
}
#endif

/**
 * @brief  Ensure the NUMA node on which we initialise (first touch) everything
 * doesn't change before engine_init allocates NUMA-local workers.
 */
void engine_pin(void) {

#ifdef HAVE_SETAFFINITY
  cpu_set_t *entry_affinity = engine_entry_affinity();

  /* Share this affinity with the threadpool, it will use this even when the
   * main thread is otherwise pinned. */
  threadpool_set_affinity_mask(entry_affinity);

  int pin;
  for (pin = 0; pin < CPU_SETSIZE && !CPU_ISSET(pin, entry_affinity); ++pin)
    ;

  cpu_set_t affinity;
  CPU_ZERO(&affinity);
  CPU_SET(pin, &affinity);
  if (sched_setaffinity(0, sizeof(affinity), &affinity) != 0) {
    error("failed to set engine's affinity.");
  }
#else
  error("SWIFT was not compiled with support for pinning.");
#endif
}

/**
 * @brief Unpins the main thread.
 */
void engine_unpin(void) {
#ifdef HAVE_SETAFFINITY
  pthread_t main_thread = pthread_self();
  cpu_set_t *entry_affinity = engine_entry_affinity();
  pthread_setaffinity_np(main_thread, sizeof(*entry_affinity), entry_affinity);
#else
  error("SWIFT was not compiled with support for pinning.");
#endif
}

/**
 * @brief Define a NUMA memory placement policy of interleave across the
 * available NUMA nodes rather than having memory in the local node, which
 * means we have a lot of memory associated with the main thread NUMA node, so
 * we don't make good use of the overall memory bandwidth between nodes.
 *
 * @param rank the MPI rank, if relevant.
 * @param verbose whether to make a report about the selected NUMA nodes.
 */
void engine_numa_policies(int rank, int verbose) {

#if defined(HAVE_LIBNUMA) && defined(_GNU_SOURCE)

  /* Get our affinity mask (on entry), that defines what NUMA nodes we should
   * use. */
  cpu_set_t *entry_affinity = engine_entry_affinity();

  /* Now convert the affinity mask into NUMA nodemask. */
  struct bitmask *nodemask = numa_allocate_nodemask();
  int nnuma = numa_num_configured_nodes();

  for (unsigned long i = 0; i < CPU_SETSIZE; i++) {

    /* If in the affinity mask we set NUMA node of CPU bit. */
    if (CPU_ISSET(i, entry_affinity)) {
      int numanode = numa_node_of_cpu(i);
      numa_bitmask_setbit(nodemask, numanode);
    }
  }

  if (verbose) {
    char report[1024];
    int len = sprintf(report, "NUMA nodes in use: [");
    for (int i = 0; i < nnuma; i++) {
      if (numa_bitmask_isbitset(nodemask, i)) {
        len += sprintf(&report[len], "%d ", i);
      } else {
        len += sprintf(&report[len], ". ");
      }
    }
    sprintf(&report[len], "]");
    pretime_message("%s", report);
  }

  /* And set. */
  set_mempolicy(MPOL_INTERLEAVE, nodemask->maskp, nodemask->size + 1);
  numa_free_nodemask(nodemask);

#endif
}

/**
 * @brief init an engine struct with the necessary properties for the
 *        simulation.
 *
 * Note do not use when restarting. Engine initialisation
 * is completed by a call to engine_config().
 *
 * @param e The #engine.
 * @param s The #space in which this #runner will run.
 * @param params The parsed parameter file.
 * @param output_options Output options for snapshots.
 * @param Ngas total number of gas particles in the simulation.
 * @param Ngparts total number of gravity particles in the simulation.
 * @param Nsinks total number of sink particles in the simulation.
 * @param Nstars total number of star particles in the simulation.
 * @param Nblackholes total number of black holes in the simulation.
 * @param Nbackground_gparts Total number of background DM particles.
 * @param Nnuparts Total number of neutrino DM particles.
 * @param policy The queuing policy to use.
 * @param verbose Is this #engine talkative ?
 * @param internal_units The system of units used internally.
 * @param physical_constants The #phys_const used for this run.
 * @param cosmo The #cosmology used for this run.
 * @param hydro The #hydro_props used for this run.
 * @param entropy_floor The #entropy_floor_properties for this run.
 * @param gravity The #gravity_props used for this run.
 * @param stars The #stars_props used for this run.
 * @param black_holes The #black_holes_props used for this run.
 * @param sinks The #sink_props used for this run.
 * @param neutrinos The #neutrino_props used for this run.
 * @param feedback The #feedback_props used for this run.
 * @param mesh The #pm_mesh used for the long-range periodic forces.
 * @param pow_data The properties and pointers for power spectrum calculation.
 * @param potential The properties of the external potential.
 * @param cooling_func The properties of the cooling function.
 * @param starform The #star_formation model of this run.
 * @param chemistry The chemistry information.
 * @param io_extra_props The properties needed for the extra i/o fields.
 * @param fof_properties The #fof_props of this run.
 * @param los_properties the #los_props of this run.
 * @param lightcone_array_properties the #lightcone_array_props of this run.
 * @param ics_metadata metadata read from the simulation ICs
 */
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
    struct ic_info *ics_metadata) {

  struct clocks_time tic, toc;
  if (engine_rank == 0) clocks_gettime(&tic);

  /* Clean-up everything */
  bzero(e, sizeof(struct engine));

  /* Store the all values in the fields of the engine. */
  e->s = s;
  e->policy = policy;
  e->step = 0;
  e->total_nr_parts = Ngas;
  e->total_nr_gparts = Ngparts;
  e->total_nr_sparts = Nstars;
  e->total_nr_sinks = Nsinks;
  e->total_nr_bparts = Nblackholes;
  e->total_nr_DM_background_gparts = Nbackground_gparts;
  e->total_nr_neutrino_gparts = Nnuparts;
  e->proxy_ind = NULL;
  e->nr_proxies = 0;
  e->ti_old = 0;
  e->ti_current = 0;
  e->ti_earliest_undrifted = 0;
  e->time_step = 0.;
  e->time_base = 0.;
  e->time_base_inv = 0.;
  e->time_begin = 0.;
  e->time_end = 0.;
  e->max_active_bin = num_time_bins;
  e->min_active_bin = 1;
  e->ti_current_subcycle = 0;
  e->max_active_bin_subcycle = num_time_bins;
  e->min_active_bin_subcycle = 1;
  e->internal_units = internal_units;
  e->output_list_snapshots = NULL;
  if (num_snapshot_triggers_part)
    parser_get_param_double_array(params, "Snapshots:recording_triggers_part",
                                  num_snapshot_triggers_part,
                                  e->snapshot_recording_triggers_desired_part);
  if (num_snapshot_triggers_spart)
    parser_get_param_double_array(params, "Snapshots:recording_triggers_spart",
                                  num_snapshot_triggers_spart,
                                  e->snapshot_recording_triggers_desired_spart);
  if (num_snapshot_triggers_bpart)
    parser_get_param_double_array(params, "Snapshots:recording_triggers_bpart",
                                  num_snapshot_triggers_bpart,
                                  e->snapshot_recording_triggers_desired_bpart);
  e->a_first_snapshot =
      parser_get_opt_param_double(params, "Snapshots:scale_factor_first", 0.1);
  e->time_first_snapshot =
      parser_get_opt_param_double(params, "Snapshots:time_first", 0.);
  e->delta_time_snapshot =
      parser_get_opt_param_double(params, "Snapshots:delta_time", -1.);
  e->ti_next_snapshot = 0;
  parser_get_param_string(params, "Snapshots:basename", e->snapshot_base_name);
  parser_get_opt_param_string(params, "Snapshots:subdir", e->snapshot_subdir,
                              engine_default_snapshot_subdir);
  parser_get_opt_param_int_array(params, "Snapshots:subsample",
                                 swift_type_count, e->snapshot_subsample);
  parser_get_opt_param_float_array(params, "Snapshots:subsample_fraction",
                                   swift_type_count,
                                   e->snapshot_subsample_fraction);
  e->snapshot_run_on_dump =
      parser_get_opt_param_int(params, "Snapshots:run_on_dump", 0);
  if (e->snapshot_run_on_dump) {
    parser_get_param_string(params, "Snapshots:dump_command",
                            e->snapshot_dump_command);
  }
  e->snapshot_compression =
      parser_get_opt_param_int(params, "Snapshots:compression", 0);
  e->snapshot_distributed =
      parser_get_opt_param_int(params, "Snapshots:distributed", 0);
  e->snapshot_lustre_OST_count =
      parser_get_opt_param_int(params, "Snapshots:lustre_OST_count", 0);
  e->snapshot_invoke_stf =
      parser_get_opt_param_int(params, "Snapshots:invoke_stf", 0);
  e->snapshot_invoke_fof =
      parser_get_opt_param_int(params, "Snapshots:invoke_fof", 0);
  e->snapshot_invoke_ps =
      parser_get_opt_param_int(params, "Snapshots:invoke_ps", 0);
  e->snapshot_use_delta_from_edge =
      parser_get_opt_param_int(params, "Snapshots:use_delta_from_edge", 0);
  if (e->snapshot_use_delta_from_edge) {
    e->snapshot_delta_from_edge =
        parser_get_param_double(params, "Snapshots:delta_from_edge");
  }
  e->dump_catalogue_when_seeding =
      parser_get_opt_param_int(params, "FOF:dump_catalogue_when_seeding", 0);
  e->snapshot_units = (struct unit_system *)malloc(sizeof(struct unit_system));
  units_init_default(e->snapshot_units, params, "Snapshots", internal_units);
  e->free_foreign_when_dumping_restart = parser_get_opt_param_int(
      params, "Scheduler:free_foreign_during_restart", 0);
  e->free_foreign_when_rebuilding = parser_get_opt_param_int(
      params, "Scheduler:free_foreign_during_rebuild", 0);
  e->snapshot_output_count = 0;
  e->stf_output_count = 0;
  e->los_output_count = 0;
  e->ps_output_count = 0;
  e->dt_min = parser_get_param_double(params, "TimeIntegration:dt_min");
  e->dt_max = parser_get_param_double(params, "TimeIntegration:dt_max");
  e->max_nr_rt_subcycles = parser_get_opt_param_int(
      params, "TimeIntegration:max_nr_rt_subcycles", /*default=*/0);
  e->dt_max_RMS_displacement = FLT_MAX;
  e->max_RMS_displacement_factor = parser_get_opt_param_double(
      params, "TimeIntegration:max_dt_RMS_factor", 0.25);
  e->max_RMS_dt_use_only_gas = parser_get_opt_param_int(
      params, "TimeIntegration:dt_RMS_use_gas_only", 0);
  e->dt_kick_grav_mesh_for_io = 0.f;
  e->a_first_statistics =
      parser_get_opt_param_double(params, "Statistics:scale_factor_first", 0.1);
  e->time_first_statistics =
      parser_get_opt_param_double(params, "Statistics:time_first", 0.);
  e->delta_time_statistics =
      parser_get_param_double(params, "Statistics:delta_time");
  e->ti_next_stats = 0;
  e->ti_next_stf = 0;
  e->ti_next_fof = 0;
  e->ti_next_ps = 0;
  e->verbose = verbose;
  e->wallclock_time = 0.f;
  e->physical_constants = physical_constants;
  e->cosmology = cosmo;
  e->hydro_properties = hydro;
  e->entropy_floor = entropy_floor;
  e->gravity_properties = gravity;
  e->stars_properties = stars;
  e->black_holes_properties = black_holes;
  e->sink_properties = sinks;
  e->neutrino_properties = neutrinos;
  e->neutrino_response = neutrino_response;
  e->mesh = mesh;
  e->power_data = pow_data;
  e->external_potential = potential;
  e->forcing_terms = forcing_terms;
  e->cooling_func = cooling_func;
  e->star_formation = starform;
  e->feedback_props = feedback;
  e->pressure_floor_props = pressure_floor;
  e->rt_props = rt;
  e->chemistry = chemistry;
  e->io_extra_props = io_extra_props;
  e->fof_properties = fof_properties;
  e->parameter_file = params;
  e->output_options = output_options;
  e->stf_this_timestep = 0;
  e->los_properties = los_properties;
  e->lightcone_array_properties = lightcone_array_properties;
  e->ics_metadata = ics_metadata;
#ifdef WITH_MPI
  e->usertime_last_step = 0.0;
  e->systime_last_step = 0.0;
  e->last_repartition = 0;
#endif
  e->total_nr_cells = 0;
  e->total_nr_tasks = 0;

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
  e->force_checks_only_all_active =
      parser_get_opt_param_int(params, "ForceChecks:only_when_all_active", 0);
  e->force_checks_only_at_snapshots =
      parser_get_opt_param_int(params, "ForceChecks:only_at_snapshots", 0);
#endif

  /* Make the space link back to the engine. */
  s->e = e;

  /* Read the run label */
  memset(e->run_name, 0, PARSER_MAX_LINE_SIZE);
  parser_get_opt_param_string(params, "MetaData:run_name", e->run_name,
                              "Untitled SWIFT simulation");
  if (strlen(e->run_name) == 0) {
    error("The run name in the parameter file cannot be an empty string.");
  }

  /* Setup the timestep if non-cosmological */
  if (!(e->policy & engine_policy_cosmology)) {
    e->time_begin =
        parser_get_param_double(params, "TimeIntegration:time_begin");
    e->time_end = parser_get_param_double(params, "TimeIntegration:time_end");
    e->time_old = e->time_begin;
    e->time = e->time_begin;

    e->time_base = (e->time_end - e->time_begin) / max_nr_timesteps;
    e->time_base_inv = 1.0 / e->time_base;
    e->ti_current = 0;
  } else {

    e->time_begin = e->cosmology->time_begin;
    e->time_end = e->cosmology->time_end;
    e->time_old = e->time_begin;
    e->time = e->time_begin;

    /* Copy the relevent information from the cosmology model */
    e->time_base = e->cosmology->time_base;
    e->time_base_inv = e->cosmology->time_base_inv;
    e->ti_current = 0;
  }

  /* Initialise VELOCIraptor output. */
  if (e->policy & engine_policy_structure_finding) {
    parser_get_param_string(params, "StructureFinding:basename",
                            e->stf_base_name);
    parser_get_opt_param_string(params, "StructureFinding:subdir_per_output",
                                e->stf_subdir_per_output,
                                engine_default_stf_subdir_per_output);
    parser_get_param_string(params, "StructureFinding:config_file_name",
                            e->stf_config_file_name);

    e->time_first_stf_output =
        parser_get_opt_param_double(params, "StructureFinding:time_first", 0.);
    e->a_first_stf_output = parser_get_opt_param_double(
        params, "StructureFinding:scale_factor_first", 0.1);
    e->delta_time_stf =
        parser_get_opt_param_double(params, "StructureFinding:delta_time", -1.);
  }

  /* Initialise line of sight output. */
  if (e->policy & engine_policy_line_of_sight) {
    e->time_first_los =
        parser_get_opt_param_double(params, "LineOfSight:time_first", 0.);
    e->a_first_los = parser_get_opt_param_double(
        params, "LineOfSight:scale_factor_first", 0.1);
    e->delta_time_los =
        parser_get_opt_param_double(params, "LineOfSight:delta_time", -1.);
  }

  /* Initialise power spectrum output. */
  if (e->policy & engine_policy_power_spectra) {
    e->time_first_ps_output =
        parser_get_opt_param_double(params, "PowerSpectrum:time_first", 0.);
    e->a_first_ps_output = parser_get_opt_param_double(
        params, "PowerSpectrum:scale_factor_first", 0.1);
    e->delta_time_ps =
        parser_get_opt_param_double(params, "PowerSpectrum:delta_time", -1.);
  }

  /* Initialise FoF calls frequency. */
  if (e->policy & engine_policy_fof) {

    e->time_first_fof_call =
        parser_get_opt_param_double(params, "FOF:time_first", 0.);
    e->a_first_fof_call =
        parser_get_opt_param_double(params, "FOF:scale_factor_first", 0.1);
    e->delta_time_fof =
        parser_get_opt_param_double(params, "FOF:delta_time", -1.);
  } else {
    if (e->snapshot_invoke_fof)
      error("Error: Must run with --fof if Snapshots::invoke_fof=1\n");
  }

  /* Initialize the star formation history structure */
  if (e->policy & engine_policy_star_formation) {
    star_formation_logger_accumulator_init(&e->sfh);
  }

  /* Initialize the neutrino mass conversion factor */
  if (Nnuparts > 0) {
    const double neutrino_volume = s->dim[0] * s->dim[1] * s->dim[2];
    e->neutrino_mass_conversion_factor =
        neutrino_mass_factor(cosmo, internal_units, physical_constants,
                             neutrino_volume, e->total_nr_neutrino_gparts);
  } else {
    e->neutrino_mass_conversion_factor = 0.f;
  }

  if (engine_rank == 0) {
    clocks_gettime(&toc);
    message("took %.3f %s.", clocks_diff(&tic, &toc), clocks_getunit());
    fflush(stdout);
  }

  /* Initialize the CSDS (already timed, not need to include it) */
#if defined(WITH_CSDS)
  if (e->policy & engine_policy_csds) {
    e->csds = (struct csds_writer *)malloc(sizeof(struct csds_writer));
    csds_init(e->csds, e, params);
  }
#endif
}

/**
 * @brief Prints the current policy of an engine
 *
 * @param e The engine to print information about
 */
void engine_print_policy(struct engine *e) {

#ifdef WITH_MPI
  if (e->nodeID == 0) {
    printf("[0000] %s engine_policy: engine policies are [ ",
           clocks_get_timesincestart());
    for (int k = 0; k <= engine_maxpolicy; k++)
      if (e->policy & (1 << k)) printf(" '%s' ", engine_policy_names[k + 1]);
    printf(" ]\n");
    fflush(stdout);
  }
#else
  printf("%s engine_policy: engine policies are [ ",
         clocks_get_timesincestart());
  for (int k = 0; k <= engine_maxpolicy; k++)
    if (e->policy & (1 << k)) printf(" '%s' ", engine_policy_names[k + 1]);
  printf(" ]\n");
  fflush(stdout);
#endif
}

/**
 * @brief Computes the maximal time-step allowed by the max RMS displacement
 * condition.
 *
 * @param e The #engine.
 */
void engine_recompute_displacement_constraint(struct engine *e) {

  const ticks tic = getticks();

  /* Get the cosmological information */
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const struct cosmology *cosmo = e->cosmology;
  const float Ocdm = cosmo->Omega_cdm;
  const float Ob = cosmo->Omega_b;
  const float H0 = cosmo->H0;
  const float a = cosmo->a;
  const float G_newton = e->physical_constants->const_newton_G;
  const float rho_crit0 = 3.f * H0 * H0 / (8.f * M_PI * G_newton);

  if (with_cosmology) {

    /* Start by reducing the minimal mass of each particle type */
    float min_mass[swift_type_count] = {e->s->min_part_mass,
                                        e->s->min_gpart_mass,
                                        FLT_MAX,
                                        e->s->min_sink_mass,
                                        e->s->min_spart_mass,
                                        e->s->min_bpart_mass,
                                        FLT_MAX};

#ifdef WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, min_mass, swift_type_count, MPI_FLOAT, MPI_MIN,
                  MPI_COMM_WORLD);
#endif

#ifdef SWIFT_DEBUG_CHECKS
    /* Check that the minimal mass collection worked */
    float min_part_mass_check = FLT_MAX;
    for (size_t i = 0; i < e->s->nr_parts; ++i) {
      if (e->s->parts[i].time_bin >= num_time_bins) continue;
      min_part_mass_check =
          min(min_part_mass_check, hydro_get_mass(&e->s->parts[i]));
    }
    if (min_part_mass_check < min_mass[swift_type_gas])
      error("Error collecting minimal mass of gas particles.");
#endif

    /* Do the same for the velocity norm sum */
    float vel_norm[swift_type_count] = {e->s->sum_part_vel_norm,
                                        e->s->sum_gpart_vel_norm,
                                        0.f,
                                        e->s->sum_sink_vel_norm,
                                        e->s->sum_spart_vel_norm,
                                        e->s->sum_spart_vel_norm,
                                        0.f};
#ifdef WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, vel_norm, swift_type_count, MPI_FLOAT, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    /* Get the counts of each particle types */
    const long long total_nr_baryons =
        e->total_nr_parts + e->total_nr_sparts + e->total_nr_bparts;
    const long long total_nr_dm_gparts =
        e->total_nr_gparts - e->total_nr_DM_background_gparts -
        e->total_nr_neutrino_gparts - total_nr_baryons;
    float count_parts[swift_type_count] = {
        (float)e->total_nr_parts,
        (float)total_nr_dm_gparts,
        (float)e->total_nr_DM_background_gparts,
        (float)e->total_nr_sinks,
        (float)e->total_nr_sparts,
        (float)e->total_nr_bparts,
        (float)e->total_nr_neutrino_gparts};

    /* Count of particles for the two species */
    const float N_dm = count_parts[1];
    const float N_b =
        count_parts[0] + count_parts[3] + count_parts[4] + count_parts[5];

    /* Peculiar motion norm for the two species */
    const float vel_norm_dm = vel_norm[1];
    const float vel_norm_b =
        vel_norm[0] + vel_norm[3] + vel_norm[4] + vel_norm[5];

    /* Mesh forces smoothing scale */
    float r_s;
    if ((e->policy & engine_policy_self_gravity) && e->s->periodic)
      r_s = e->mesh->r_s;
    else
      r_s = FLT_MAX;

    float dt_dm = FLT_MAX, dt_b = FLT_MAX;

    /* DM case */
    if (N_dm > 0.f) {

      /* Minimal mass for the DM */
      const float min_mass_dm = min_mass[1];

      /* Inter-particle sepration for the DM */
      const float d_dm = cbrtf(min_mass_dm / (Ocdm * rho_crit0));

      /* RMS peculiar motion for the DM */
      const float rms_vel_dm = vel_norm_dm / N_dm;

      /* Time-step based on maximum displacement */
      dt_dm = a * a * min(r_s, d_dm) / sqrtf(rms_vel_dm);
    }

    /* Baryon case */
    if (N_b > 0.f) {

      /* Minimal mass for the bayons */
      float min_mass_b;
      if (e->max_RMS_dt_use_only_gas)
        min_mass_b = min_mass[0];
      else
        min_mass_b = min4(min_mass[0], min_mass[3], min_mass[4], min_mass[5]);

      /* Inter-particle sepration for the baryons */
      const float d_b = cbrtf(min_mass_b / (Ob * rho_crit0));

      /* RMS peculiar motion for the baryons */
      const float rms_vel_b = vel_norm_b / N_b;

      /* Time-step based on maximum displacement */
      dt_b = a * a * min(r_s, d_b) / sqrtf(rms_vel_b);
    }

    /* Use the minimum */
    const float dt = min(dt_dm, dt_b);

    /* Apply the dimensionless factor */
    e->dt_max_RMS_displacement = dt * e->max_RMS_displacement_factor;

    if (e->verbose)
      message("max_dt_RMS_displacement = %e", e->dt_max_RMS_displacement);
  }

  /* Now, update the mesh time-step */

  /* Store the previous time-step size */
  e->mesh->ti_end_mesh_last = e->mesh->ti_end_mesh_next;
  e->mesh->ti_beg_mesh_last = e->mesh->ti_beg_mesh_next;
  const integertime_t old_dti =
      e->mesh->ti_end_mesh_last - e->mesh->ti_beg_mesh_last;

#ifdef SWIFT_DEBUG_CHECKS
  if (e->step > 1 && e->mesh->ti_end_mesh_last != e->ti_current)
    error("Weird time integration issue");
#endif

  /* What is the allowed time-step size
   * Note: The cosmology factor is 1 in non-cosmo runs */
  double dt_mesh = e->dt_max_RMS_displacement * e->cosmology->time_step_factor;
  dt_mesh = min(dt_mesh, e->dt_max);

  /* Convert to integer time */
  integertime_t new_dti = (integertime_t)(dt_mesh * e->time_base_inv);

  /* Find the max integer time-step on the timeline below new_dti */
  integertime_t dti_timeline = max_nr_timesteps;
  while (new_dti < dti_timeline) dti_timeline /= ((integertime_t)2);
  new_dti = dti_timeline;

  /* Make sure we are allowed to increase the timestep size */
  const integertime_t current_dti = e->step > 0 ? old_dti : max_nr_timesteps;
  if (new_dti > current_dti) {
    if ((max_nr_timesteps - e->ti_current) % new_dti > 0) {
      new_dti = current_dti;
    }
  }

  e->mesh->ti_beg_mesh_next = e->ti_current;
  e->mesh->ti_end_mesh_next = e->ti_current + new_dti;

  const timebin_t bin = get_time_bin(new_dti);

  if (e->verbose && new_dti != old_dti)
    message("Mesh time-step changed to %e (time-bin %d)",
            get_timestep(bin, e->time_base), bin);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Frees up the memory allocated for this #engine
 *
 * @param e The #engine to clean.
 * @param fof Was this a stand-alone FOF run?
 * @param restart Was this a run that was restarted from check-point files?
 */
void engine_clean(struct engine *e, const int fof, const int restart) {
  /* Start by telling the runners to stop. */
  e->step_props = engine_step_prop_done;
  swift_barrier_wait(&e->run_barrier);

  /* Wait for each runner to come home. */
  for (int k = 0; k < e->nr_threads; k++) {
    if (pthread_join(e->runners[k].thread, /*retval=*/NULL) != 0)
      error("Failed to join runner %i.", k);
#ifdef WITH_VECTORIZATION
    cache_clean(&e->runners[k].ci_cache);
    cache_clean(&e->runners[k].cj_cache);
#endif
    gravity_cache_clean(&e->runners[k].ci_gravity_cache);
    gravity_cache_clean(&e->runners[k].cj_gravity_cache);
  }
  swift_free("runners", e->runners);
  free(e->snapshot_units);

  output_list_clean(&e->output_list_snapshots);
  output_list_clean(&e->output_list_stats);
  output_list_clean(&e->output_list_stf);
  output_list_clean(&e->output_list_los);
  output_list_clean(&e->output_list_ps);

  output_options_clean(e->output_options);

  ic_info_clean(e->ics_metadata);

  swift_free("links", e->links);
#if defined(WITH_CSDS)
  if (e->policy & engine_policy_csds) {
    csds_free(e->csds);
    free(e->csds);
  }
#endif
  scheduler_clean(&e->sched);
  space_clean(e->s);
  threadpool_clean(&e->threadpool);
#if defined(WITH_MPI)
  for (int i = 0; i < e->nr_proxies; ++i) {
    proxy_clean(&e->proxies[i]);
  }
  free(e->proxy_ind);
  free(e->proxies);

  /* Free types */
  part_free_mpi_types();
  multipole_free_mpi_types();
  stats_free_mpi_type();
  proxy_free_mpi_type();
  task_free_mpi_comms();
  if (!fof) mpicollect_free_MPI_type();
#endif

  /* Close files */
  if (!fof && e->nodeID == 0) {
    fclose(e->file_timesteps);
    fclose(e->file_stats);

    if (e->policy & engine_policy_star_formation) {
      fclose(e->sfh_logger);
    }

#ifndef RT_NONE
    fclose(e->file_rt_subcycles);
#endif
  }

  /* If the run was restarted, we should also free the memory allocated
     in engine_struct_restore() */
  if (restart) {
    free((void *)e->parameter_file);
    free((void *)e->output_options);
    free((void *)e->external_potential);
    free((void *)e->forcing_terms);
    free((void *)e->black_holes_properties);
    free((void *)e->pressure_floor_props);
    free((void *)e->rt_props);
    free((void *)e->sink_properties);
    free((void *)e->stars_properties);
    free((void *)e->gravity_properties);
    free((void *)e->neutrino_properties);
    free((void *)e->hydro_properties);
    free((void *)e->physical_constants);
    free((void *)e->internal_units);
    free((void *)e->cosmology);
    free((void *)e->mesh);
    free((void *)e->power_data);
    free((void *)e->chemistry);
    free((void *)e->entropy_floor);
    free((void *)e->cooling_func);
    free((void *)e->star_formation);
    free((void *)e->feedback_props);
    free((void *)e->io_extra_props);
#ifdef WITH_FOF
    free((void *)e->fof_properties);
#endif
    free((void *)e->los_properties);
    free((void *)e->lightcone_array_properties);
    free((void *)e->ics_metadata);
#ifdef WITH_MPI
    free((void *)e->reparttype);
#endif
    if (e->output_list_snapshots) free((void *)e->output_list_snapshots);
    if (e->output_list_stats) free((void *)e->output_list_stats);
    if (e->output_list_stf) free((void *)e->output_list_stf);
    if (e->output_list_los) free((void *)e->output_list_los);
    if (e->output_list_ps) free((void *)e->output_list_ps);
#ifdef WITH_CSDS
    if (e->policy & engine_policy_csds) free((void *)e->csds);
#endif
    free(e->s);
  }
}

/**
 * @brief Write the engine struct and its contents to the given FILE as a
 * stream of bytes.
 *
 * @param e the engine
 * @param stream the file stream
 */
void engine_struct_dump(struct engine *e, FILE *stream) {

  /* Dump the engine. Save the current tasks_per_cell estimate. */
  e->restart_max_tasks = engine_estimate_nr_tasks(e);
  restart_write_blocks(e, sizeof(struct engine), 1, stream, "engine",
                       "engine struct");

  /* And all the engine pointed data, these use their own dump functions. */
  space_struct_dump(e->s, stream);
  units_struct_dump(e->internal_units, stream);
  units_struct_dump(e->snapshot_units, stream);
  cosmology_struct_dump(e->cosmology, stream);

#ifdef WITH_MPI
  /* Save the partition for restoration. */
  partition_store_celllist(e->s, e->reparttype);
  partition_struct_dump(e->reparttype, stream);
#endif

  phys_const_struct_dump(e->physical_constants, stream);
  hydro_props_struct_dump(e->hydro_properties, stream);
  entropy_floor_struct_dump(e->entropy_floor, stream);
  gravity_props_struct_dump(e->gravity_properties, stream);
  stars_props_struct_dump(e->stars_properties, stream);
  pm_mesh_struct_dump(e->mesh, stream);
  power_spectrum_struct_dump(e->power_data, stream);
  potential_struct_dump(e->external_potential, stream);
  forcing_terms_struct_dump(e->forcing_terms, stream);
  cooling_struct_dump(e->cooling_func, stream);
  starformation_struct_dump(e->star_formation, stream);
  feedback_struct_dump(e->feedback_props, stream);
  pressure_floor_struct_dump(e->pressure_floor_props, stream);
  rt_struct_dump(e->rt_props, stream);
  black_holes_struct_dump(e->black_holes_properties, stream);
  sink_struct_dump(e->sink_properties, stream);
  neutrino_struct_dump(e->neutrino_properties, stream);
  neutrino_response_struct_dump(e->neutrino_response, stream);
  chemistry_struct_dump(e->chemistry, stream);
  extra_io_struct_dump(e->io_extra_props, stream);
#ifdef WITH_FOF
  fof_struct_dump(e->fof_properties, stream);
#endif
  los_struct_dump(e->los_properties, stream);
  lightcone_array_struct_dump(e->lightcone_array_properties, stream);
  ic_info_struct_dump(e->ics_metadata, stream);
  parser_struct_dump(e->parameter_file, stream);
  output_options_struct_dump(e->output_options, stream);

#ifdef WITH_CSDS
  if (e->policy & engine_policy_csds) {
    csds_struct_dump(e->csds, stream);
  }
#endif
}

/**
 * @brief Re-create an engine struct and its contents from the given FILE
 *        stream.
 *
 * @param e the engine
 * @param stream the file stream
 */
void engine_struct_restore(struct engine *e, FILE *stream) {

  /* Read the engine. */
  restart_read_blocks(e, sizeof(struct engine), 1, stream, NULL,
                      "engine struct");

  /* Re-initializations as necessary for our struct and its members. */
  e->sched.tasks = NULL;
  e->sched.tasks_ind = NULL;
  e->sched.tid_active = NULL;
  e->sched.size = 0;

  /* Now for the other pointers, these use their own restore functions. */
  /* Note all this memory leaks, but is used once. */
  struct space *s = (struct space *)malloc(sizeof(struct space));
  space_struct_restore(s, stream);
  e->s = s;
  s->e = e;

  struct unit_system *internal_us =
      (struct unit_system *)malloc(sizeof(struct unit_system));
  units_struct_restore(internal_us, stream);
  e->internal_units = internal_us;

  struct unit_system *snap_us =
      (struct unit_system *)malloc(sizeof(struct unit_system));
  units_struct_restore(snap_us, stream);
  e->snapshot_units = snap_us;

  struct cosmology *cosmo =
      (struct cosmology *)malloc(sizeof(struct cosmology));
  cosmology_struct_restore(e->policy & engine_policy_cosmology, cosmo, stream);
  e->cosmology = cosmo;

#ifdef WITH_MPI
  struct repartition *reparttype =
      (struct repartition *)malloc(sizeof(struct repartition));
  partition_struct_restore(reparttype, stream);
  e->reparttype = reparttype;
#endif

  struct phys_const *physical_constants =
      (struct phys_const *)malloc(sizeof(struct phys_const));
  phys_const_struct_restore(physical_constants, stream);
  e->physical_constants = physical_constants;

  struct hydro_props *hydro_properties =
      (struct hydro_props *)malloc(sizeof(struct hydro_props));
  hydro_props_struct_restore(hydro_properties, stream);
  e->hydro_properties = hydro_properties;

  struct entropy_floor_properties *entropy_floor =
      (struct entropy_floor_properties *)malloc(
          sizeof(struct entropy_floor_properties));
  entropy_floor_struct_restore(entropy_floor, stream);
  e->entropy_floor = entropy_floor;

  struct gravity_props *gravity_properties =
      (struct gravity_props *)malloc(sizeof(struct gravity_props));
  gravity_props_struct_restore(gravity_properties, stream);
  e->gravity_properties = gravity_properties;

  struct stars_props *stars_properties =
      (struct stars_props *)malloc(sizeof(struct stars_props));
  stars_props_struct_restore(stars_properties, stream);
  e->stars_properties = stars_properties;

  struct pm_mesh *mesh = (struct pm_mesh *)malloc(sizeof(struct pm_mesh));
  pm_mesh_struct_restore(mesh, stream);
  e->mesh = mesh;

  struct power_spectrum_data *pow_data =
      (struct power_spectrum_data *)malloc(sizeof(struct power_spectrum_data));
  power_spectrum_struct_restore(pow_data, stream);
  e->power_data = pow_data;

  struct external_potential *external_potential =
      (struct external_potential *)malloc(sizeof(struct external_potential));
  potential_struct_restore(external_potential, stream);
  e->external_potential = external_potential;

  struct forcing_terms *forcing_terms =
      (struct forcing_terms *)malloc(sizeof(struct forcing_terms));
  forcing_terms_struct_restore(forcing_terms, stream);
  e->forcing_terms = forcing_terms;

  struct cooling_function_data *cooling_func =
      (struct cooling_function_data *)malloc(
          sizeof(struct cooling_function_data));
  cooling_struct_restore(cooling_func, stream, e->cosmology);
  e->cooling_func = cooling_func;

  struct star_formation *star_formation =
      (struct star_formation *)malloc(sizeof(struct star_formation));
  starformation_struct_restore(star_formation, stream);
  e->star_formation = star_formation;

  struct feedback_props *feedback_properties =
      (struct feedback_props *)malloc(sizeof(struct feedback_props));
  feedback_struct_restore(feedback_properties, stream);
  e->feedback_props = feedback_properties;

  struct pressure_floor_props *pressure_floor_properties =
      (struct pressure_floor_props *)malloc(
          sizeof(struct pressure_floor_props));
  pressure_floor_struct_restore(pressure_floor_properties, stream);
  e->pressure_floor_props = pressure_floor_properties;

  struct rt_props *rt_properties =
      (struct rt_props *)malloc(sizeof(struct rt_props));
  rt_struct_restore(rt_properties, stream, e->physical_constants,
                    e->internal_units);
  e->rt_props = rt_properties;

  struct black_holes_props *black_holes_properties =
      (struct black_holes_props *)malloc(sizeof(struct black_holes_props));
  black_holes_struct_restore(black_holes_properties, stream);
  e->black_holes_properties = black_holes_properties;

  struct sink_props *sink_properties =
      (struct sink_props *)malloc(sizeof(struct sink_props));
  sink_struct_restore(sink_properties, stream);
  e->sink_properties = sink_properties;

  struct neutrino_props *neutrino_properties =
      (struct neutrino_props *)malloc(sizeof(struct neutrino_props));
  neutrino_struct_restore(neutrino_properties, stream);
  e->neutrino_properties = neutrino_properties;

  struct neutrino_response *neutrino_response =
      (struct neutrino_response *)malloc(sizeof(struct neutrino_response));
  neutrino_response_struct_restore(neutrino_response, stream);
  e->neutrino_response = neutrino_response;

  struct chemistry_global_data *chemistry =
      (struct chemistry_global_data *)malloc(
          sizeof(struct chemistry_global_data));
  chemistry_struct_restore(chemistry, stream);
  e->chemistry = chemistry;

  struct extra_io_properties *extra_io_props =
      (struct extra_io_properties *)malloc(sizeof(struct extra_io_properties));
  extra_io_struct_restore(extra_io_props, stream);
  e->io_extra_props = extra_io_props;

#ifdef WITH_FOF
  struct fof_props *fof_props =
      (struct fof_props *)malloc(sizeof(struct fof_props));
  fof_struct_restore(fof_props, stream);
  e->fof_properties = fof_props;
#endif

  struct los_props *los_properties =
      (struct los_props *)malloc(sizeof(struct los_props));
  los_struct_restore(los_properties, stream);
  e->los_properties = los_properties;

  struct lightcone_array_props *lightcone_array_properties =
      (struct lightcone_array_props *)malloc(
          sizeof(struct lightcone_array_props));
  lightcone_array_struct_restore(lightcone_array_properties, stream);
  e->lightcone_array_properties = lightcone_array_properties;

  struct ic_info *ics_metadata =
      (struct ic_info *)malloc(sizeof(struct ic_info));
  ic_info_struct_restore(ics_metadata, stream);
  e->ics_metadata = ics_metadata;

  struct swift_params *parameter_file =
      (struct swift_params *)malloc(sizeof(struct swift_params));
  parser_struct_restore(parameter_file, stream);
  e->parameter_file = parameter_file;

  struct output_options *output_options =
      (struct output_options *)malloc(sizeof(struct output_options));
  output_options_struct_restore(output_options, stream);
  e->output_options = output_options;

#ifdef WITH_CSDS
  if (e->policy & engine_policy_csds) {
    struct csds_writer *log =
        (struct csds_writer *)malloc(sizeof(struct csds_writer));
    csds_struct_restore(log, stream);
    e->csds = log;
  }
#endif

#ifdef EOS_PLANETARY
  eos_init(&eos, e->physical_constants, e->snapshot_units, e->parameter_file);
#endif

  /* Want to force a rebuild before using this engine. Wait to repartition.*/
  e->forcerebuild = 1;
  e->forcerepart = 0;
}
