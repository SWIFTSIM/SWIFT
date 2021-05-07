/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "csds_io.h"
#include "distributed_io.h"
#include "kick.h"
#include "lightcone.h"
#include "lightcone_array.h"
#include "line_of_sight.h"
#include "parallel_io.h"
#include "serial_io.h"
#include "single_io.h"

#include <stdio.h>

/**
 * @brief dump restart files if it is time to do so and dumps are enabled.
 *
 * @param e the engine.
 * @param drifted_all true if a drift_all has just been performed.
 * @param force force a dump, if dumping is enabled.
 */
void engine_dump_restarts(struct engine *e, int drifted_all, int force) {

  if (e->restart_dump) {
    ticks tic = getticks();

    /* Dump when the time has arrived, or we are told to. */
    int dump = ((tic > e->restart_next) || force);

#ifdef WITH_MPI
    /* Synchronize this action from rank 0 (ticks may differ between
     * machines). */
    MPI_Bcast(&dump, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if (dump) {

      if (e->nodeID == 0) message("Writing restart files");

      /* Clean out the previous saved files, if found. Do this now as we are
       * MPI synchronized. */
      restart_remove_previous(e->restart_file);

      /* Drift all particles first (may have just been done). */
      if (!drifted_all) engine_drift_all(e, /*drift_mpole=*/1);

#ifdef WITH_LIGHTCONE 
      /* Flush lightcone buffers before dumping restarts */
      lightcone_array_flush(e->lightcone_array_properties,
                            e->cosmology, e->internal_units, e->snapshot_units,
                            /*flush_map_updates=*/1, /*flush_particles=*/1,
                            /*end_file=*/1, /*dump_all_shells=*/0);
#ifdef WITH_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

      restart_write(e, e->restart_file);

#ifdef WITH_MPI
      /* Make sure all ranks finished writing to avoid having incomplete
       * sets of restart files should the code crash before all the ranks
       * are done */
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      if (e->verbose)
        message("Dumping restart files took %.3f %s",
                clocks_from_ticks(getticks() - tic), clocks_getunit());

      /* Time after which next dump will occur. */
      e->restart_next += e->restart_dt;

      /* Flag that we dumped the restarts */
      e->step_props |= engine_step_prop_restarts;
    }
  }
}

/**
 * @brief Writes a snapshot with the current state of the engine
 *
 * @param e The #engine.
 */
void engine_dump_snapshot(struct engine *e) {

  struct clocks_time time1, time2;
  clocks_gettime(&time1);

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that all cells have been drifted to the current time.
   * That can include cells that have not
   * previously been active on this rank. */
  space_check_drift_point(e->s, e->ti_current, /* check_mpole=*/0);

  /* Be verbose about this */
  if (e->nodeID == 0) {
    if (e->policy & engine_policy_cosmology)
      message("Dumping snapshot at a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Dumping snapshot at t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }
#else
  if (e->verbose) {
    if (e->policy & engine_policy_cosmology)
      message("Dumping snapshot at a=%e",
              exp(e->ti_current * e->time_base) * e->cosmology->a_begin);
    else
      message("Dumping snapshot at t=%e",
              e->ti_current * e->time_base + e->time_begin);
  }
#endif

#ifdef DEBUG_INTERACTIONS_STARS
  engine_collect_stars_counter(e);
#endif

  /* Get time-step since the last mesh kick */
  if ((e->policy & engine_policy_self_gravity) && e->s->periodic) {
    const int with_cosmology = e->policy & engine_policy_cosmology;

    e->dt_kick_grav_mesh_for_io =
        kick_get_grav_kick_dt(e->mesh->ti_beg_mesh_next, e->ti_current,
                              e->time_base, with_cosmology, e->cosmology) -
        kick_get_grav_kick_dt(
            e->mesh->ti_beg_mesh_next,
            (e->mesh->ti_beg_mesh_next + e->mesh->ti_end_mesh_next) / 2,
            e->time_base, with_cosmology, e->cosmology);
  }

/* Dump (depending on the chosen strategy) ... */
#if defined(HAVE_HDF5)
#if defined(WITH_MPI)

  if (e->snapshot_distributed) {

    write_output_distributed(e, e->internal_units, e->snapshot_units, e->nodeID,
                             e->nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
  } else {

#if defined(HAVE_PARALLEL_HDF5)
    write_output_parallel(e, e->internal_units, e->snapshot_units, e->nodeID,
                          e->nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
    write_output_serial(e, e->internal_units, e->snapshot_units, e->nodeID,
                        e->nr_nodes, MPI_COMM_WORLD, MPI_INFO_NULL);
#endif
  }
#else
  write_output_single(e, e->internal_units, e->snapshot_units);
#endif
#endif

  /* Flag that we dumped a snapshot */
  e->step_props |= engine_step_prop_snapshot;

  clocks_gettime(&time2);
  if (e->verbose)
    message("writing particle properties took %.3f %s.",
            (float)clocks_diff(&time1, &time2), clocks_getunit());

  /* Run the post-dump command if required */
  engine_run_on_dump(e);
}

/**
 * @brief Runs the snapshot_dump_command if relevant. Note that we
 *        perform no error checking on this command, and assume
 *        it works fine.
 *
 * @param e The #engine.
 */
void engine_run_on_dump(struct engine *e) {
  if (e->snapshot_run_on_dump) {
    /* Generate a string containing (optionally) the snapshot number.
     * Note that -1 is used because snapshot_output_count was just
     * increased when the write_output_* functions are called. */
    const int buf_size = PARSER_MAX_LINE_SIZE * 3;
    char dump_command_buf[buf_size];
    snprintf(dump_command_buf, buf_size, "%s %s %04d", e->snapshot_dump_command,
             e->snapshot_base_name, e->snapshot_output_count - 1);

    /* Let's trust the user's command... */
    const int result = system(dump_command_buf);
    if (result != 0) {
      message("Snapshot dump command returned error code %d", result);
    }
  }
}

/**
 * @brief Check whether any kind of i/o has to be performed during this
 * step.
 *
 * This includes snapshots, stats and halo finder. We also handle the case
 * of multiple outputs between two steps.
 *
 * @param e The #engine.
 */
void engine_check_for_dumps(struct engine *e) {
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int with_stf = (e->policy & engine_policy_structure_finding);
  const int with_los = (e->policy & engine_policy_line_of_sight);
  const int with_fof = (e->policy & engine_policy_fof);

  /* What kind of output are we getting? */
  enum output_type {
    output_none,
    output_snapshot,
    output_statistics,
    output_stf,
    output_los,
  };

  /* What kind of output do we want? And at which time ?
   * Find the earliest output (amongst all kinds) that takes place
   * before the next time-step */
  enum output_type type = output_none;
  integertime_t ti_output = max_nr_timesteps;
  e->stf_this_timestep = 0;

  /* Save some statistics ? */
  if (e->ti_end_min > e->ti_next_stats && e->ti_next_stats > 0) {
    if (e->ti_next_stats < ti_output) {
      ti_output = e->ti_next_stats;
      type = output_statistics;
    }
  }

  /* Do we want a snapshot? */
  if (e->ti_end_min > e->ti_next_snapshot && e->ti_next_snapshot > 0) {
    if (e->ti_next_snapshot < ti_output) {
      ti_output = e->ti_next_snapshot;
      type = output_snapshot;
    }
  }

  /* Do we want to perform structure finding? */
  if (with_stf) {
    if (e->ti_end_min > e->ti_next_stf && e->ti_next_stf > 0) {
      if (e->ti_next_stf < ti_output) {
        ti_output = e->ti_next_stf;
        type = output_stf;
      }
    }
  }

  /* Do we want to write a line of sight file? */
  if (with_los) {
    if (e->ti_end_min > e->ti_next_los && e->ti_next_los > 0) {
      if (e->ti_next_los < ti_output) {
        ti_output = e->ti_next_los;
        type = output_los;
      }
    }
  }

  /* Store information before attempting extra dump-related drifts */
  const integertime_t ti_current = e->ti_current;
  const timebin_t max_active_bin = e->max_active_bin;
  const double time = e->time;

  while (type != output_none) {

    /* Let's fake that we are at the dump time */
    e->ti_current = ti_output;
    e->max_active_bin = 0;
    if (with_cosmology) {
      cosmology_update(e->cosmology, e->physical_constants, e->ti_current);
      e->time = e->cosmology->time;
    } else {
      e->time = ti_output * e->time_base + e->time_begin;
    }

    /* Drift everyone */
    engine_drift_all(e, /*drift_mpole=*/0);

    /* Write some form of output */
    switch (type) {

      case output_snapshot:

#ifdef SWIFT_GRAVITY_FORCE_CHECKS
        /* Indicate we are allowed to do a brute force calculation now */
        e->force_checks_snapshot_flag = 1;
#endif

        /* Free the mesh memory to get some breathing space */
        if ((e->policy & engine_policy_self_gravity) && e->s->periodic)
          pm_mesh_free(e->mesh);

        /* Do we want FoF group IDs in the snapshot? */
        if (with_fof && e->snapshot_invoke_fof) {
          engine_fof(e, /*dump_results=*/1, /*dump_debug=*/0,
                     /*seed_black_holes=*/0);
        }

        /* Free the foreign particles to get more breathing space.
         * This cannot be done before FOF as comms are used in there. */
#ifdef WITH_MPI
        space_free_foreign_parts(e->s, /*clear_cell_pointers=*/1);
#endif

        /* Do we want a corresponding VELOCIraptor output? */
        if (with_stf && e->snapshot_invoke_stf && !e->stf_this_timestep) {

#ifdef HAVE_VELOCIRAPTOR
          velociraptor_invoke(e, /*linked_with_snap=*/1);
          e->step_props |= engine_step_prop_stf;
#else
          error(
              "Asking for a VELOCIraptor output but SWIFT was compiled without "
              "the interface!");
#endif
        }

        /* Dump... */
        engine_dump_snapshot(e);

        /* Free the memory allocated for VELOCIraptor i/o. */
        if (with_stf && e->snapshot_invoke_stf && e->s->gpart_group_data) {
#ifdef HAVE_VELOCIRAPTOR
          swift_free("gpart_group_data", e->s->gpart_group_data);
          e->s->gpart_group_data = NULL;
#endif
        }

        /* Reallocate freed memory */
        if ((e->policy & engine_policy_self_gravity) && e->s->periodic)
          pm_mesh_allocate(e->mesh);
#ifdef WITH_MPI
        engine_allocate_foreign_particles(e);
#endif

        /* ... and find the next output time */
        engine_compute_next_snapshot_time(e);
        break;

      case output_statistics:

        /* Dump */
        engine_print_stats(e);

        /* and move on */
        engine_compute_next_statistics_time(e);

        break;

      case output_stf:

        /* Free the mesh memory to get some breathing space */
        if ((e->policy & engine_policy_self_gravity) && e->s->periodic)
          pm_mesh_free(e->mesh);
#ifdef WITH_MPI
        space_free_foreign_parts(e->s, /*clear_cell_pointers=*/1);
#endif

#ifdef HAVE_VELOCIRAPTOR
        /* Unleash the raptor! */
        if (!e->stf_this_timestep) {
          velociraptor_invoke(e, /*linked_with_snap=*/0);
          e->step_props |= engine_step_prop_stf;
        }

        /* Reallocate freed memory */
        if ((e->policy & engine_policy_self_gravity) && e->s->periodic)
          pm_mesh_allocate(e->mesh);
#ifdef WITH_MPI
        engine_allocate_foreign_particles(e);
#endif

        /* ... and find the next output time */
        engine_compute_next_stf_time(e);
#else
        error(
            "Asking for a VELOCIraptor output but SWIFT was compiled without "
            "the interface!");
#endif
        break;

      case output_los:

        /* Compute the LoS */
        do_line_of_sight(e);

        /* Move on */
        engine_compute_next_los_time(e);

        break;

      default:
        error("Invalid dump type");
    }

    /* We need to see whether whether we are in the pathological case
     * where there can be another dump before the next step. */

    type = output_none;
    ti_output = max_nr_timesteps;

    /* Save some statistics ? */
    if (e->ti_end_min > e->ti_next_stats && e->ti_next_stats > 0) {
      if (e->ti_next_stats < ti_output) {
        ti_output = e->ti_next_stats;
        type = output_statistics;
      }
    }

    /* Do we want a snapshot? */
    if (e->ti_end_min > e->ti_next_snapshot && e->ti_next_snapshot > 0) {
      if (e->ti_next_snapshot < ti_output) {
        ti_output = e->ti_next_snapshot;
        type = output_snapshot;
      }
    }

    /* Do we want to perform structure finding? */
    if (with_stf) {
      if (e->ti_end_min > e->ti_next_stf && e->ti_next_stf > 0) {
        if (e->ti_next_stf < ti_output) {
          ti_output = e->ti_next_stf;
          type = output_stf;
        }
      }
    }

    /* Do line of sight ? */
    if (with_los) {
      if (e->ti_end_min > e->ti_next_los && e->ti_next_los > 0) {
        if (e->ti_next_los < ti_output) {
          ti_output = e->ti_next_los;
          type = output_los;
        }
      }
    }

  } /* While loop over output types */

  /* Restore the information we stored */
  e->ti_current = ti_current;
  if (e->policy & engine_policy_cosmology)
    cosmology_update(e->cosmology, e->physical_constants, e->ti_current);
  e->max_active_bin = max_active_bin;
  e->time = time;
}

/**
 * @brief Computes the next time (on the time line) for a dump
 *
 * @param e The #engine.
 */
void engine_compute_next_snapshot_time(struct engine *e) {

  /* Do output_list file case */
  if (e->output_list_snapshots) {
    output_list_read_next_time(e->output_list_snapshots, e, "snapshots",
                               &e->ti_next_snapshot);
    return;
  }

  /* Find upper-bound on last output */
  double time_end;
  if (e->policy & engine_policy_cosmology)
    time_end = e->cosmology->a_end * e->delta_time_snapshot;
  else
    time_end = e->time_end + e->delta_time_snapshot;

  /* Find next snasphot above current time */
  double time;
  if (e->policy & engine_policy_cosmology)
    time = e->a_first_snapshot;
  else
    time = e->time_first_snapshot;

  int found_snapshot_time = 0;
  while (time < time_end) {

    /* Output time on the integer timeline */
    if (e->policy & engine_policy_cosmology)
      e->ti_next_snapshot = log(time / e->cosmology->a_begin) / e->time_base;
    else
      e->ti_next_snapshot = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (e->ti_next_snapshot > e->ti_current) {
      found_snapshot_time = 1;
      break;
    }

    if (e->policy & engine_policy_cosmology)
      time *= e->delta_time_snapshot;
    else
      time += e->delta_time_snapshot;
  }

  /* Deal with last snapshot */
  if (!found_snapshot_time) {
    e->ti_next_snapshot = -1;
    if (e->verbose) message("No further output time.");
  } else {

    /* Be nice, talk... */
    if (e->policy & engine_policy_cosmology) {
      const double next_snapshot_time =
          exp(e->ti_next_snapshot * e->time_base) * e->cosmology->a_begin;
      if (e->verbose)
        message("Next snapshot time set to a=%e.", next_snapshot_time);
    } else {
      const double next_snapshot_time =
          e->ti_next_snapshot * e->time_base + e->time_begin;
      if (e->verbose)
        message("Next snapshot time set to t=%e.", next_snapshot_time);
    }
  }
}

/**
 * @brief Computes the next time (on the time line) for a statistics dump
 *
 * @param e The #engine.
 */
void engine_compute_next_statistics_time(struct engine *e) {
  /* Do output_list file case */
  if (e->output_list_stats) {
    output_list_read_next_time(e->output_list_stats, e, "stats",
                               &e->ti_next_stats);
    return;
  }

  /* Find upper-bound on last output */
  double time_end;
  if (e->policy & engine_policy_cosmology)
    time_end = e->cosmology->a_end * e->delta_time_statistics;
  else
    time_end = e->time_end + e->delta_time_statistics;

  /* Find next snasphot above current time */
  double time;
  if (e->policy & engine_policy_cosmology)
    time = e->a_first_statistics;
  else
    time = e->time_first_statistics;

  int found_stats_time = 0;
  while (time < time_end) {

    /* Output time on the integer timeline */
    if (e->policy & engine_policy_cosmology)
      e->ti_next_stats = log(time / e->cosmology->a_begin) / e->time_base;
    else
      e->ti_next_stats = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (e->ti_next_stats > e->ti_current) {
      found_stats_time = 1;
      break;
    }

    if (e->policy & engine_policy_cosmology)
      time *= e->delta_time_statistics;
    else
      time += e->delta_time_statistics;
  }

  /* Deal with last statistics */
  if (!found_stats_time) {
    e->ti_next_stats = -1;
    if (e->verbose) message("No further output time.");
  } else {

    /* Be nice, talk... */
    if (e->policy & engine_policy_cosmology) {
      const double next_statistics_time =
          exp(e->ti_next_stats * e->time_base) * e->cosmology->a_begin;
      if (e->verbose)
        message("Next output time for stats set to a=%e.",
                next_statistics_time);
    } else {
      const double next_statistics_time =
          e->ti_next_stats * e->time_base + e->time_begin;
      if (e->verbose)
        message("Next output time for stats set to t=%e.",
                next_statistics_time);
    }
  }
}

/**
 * @brief Computes the next time (on the time line) for a line of sight dump
 *
 * @param e The #engine.
 */
void engine_compute_next_los_time(struct engine *e) {
  /* Do output_list file case */
  if (e->output_list_los) {
    output_list_read_next_time(e->output_list_los, e, "line of sights",
                               &e->ti_next_los);
    return;
  }

  /* Find upper-bound on last output */
  double time_end;
  if (e->policy & engine_policy_cosmology)
    time_end = e->cosmology->a_end * e->delta_time_los;
  else
    time_end = e->time_end + e->delta_time_los;

  /* Find next los above current time */
  double time;
  if (e->policy & engine_policy_cosmology)
    time = e->a_first_los;
  else
    time = e->time_first_los;

  int found_los_time = 0;
  while (time < time_end) {

    /* Output time on the integer timeline */
    if (e->policy & engine_policy_cosmology)
      e->ti_next_los = log(time / e->cosmology->a_begin) / e->time_base;
    else
      e->ti_next_los = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (e->ti_next_los > e->ti_current) {
      found_los_time = 1;
      break;
    }

    if (e->policy & engine_policy_cosmology)
      time *= e->delta_time_los;
    else
      time += e->delta_time_los;
  }

  /* Deal with last line of sight */
  if (!found_los_time) {
    e->ti_next_los = -1;
    if (e->verbose) message("No further LOS output time.");
  } else {

    /* Be nice, talk... */
    if (e->policy & engine_policy_cosmology) {
      const double next_los_time =
          exp(e->ti_next_los * e->time_base) * e->cosmology->a_begin;
      if (e->verbose)
        message("Next output time for line of sight set to a=%e.",
                next_los_time);
    } else {
      const double next_los_time =
          e->ti_next_los * e->time_base + e->time_begin;
      if (e->verbose)
        message("Next output time for line of sight set to t=%e.",
                next_los_time);
    }
  }
}

/**
 * @brief Computes the next time (on the time line) for structure finding
 *
 * @param e The #engine.
 */
void engine_compute_next_stf_time(struct engine *e) {
  /* Do output_list file case */
  if (e->output_list_stf) {
    output_list_read_next_time(e->output_list_stf, e, "stf", &e->ti_next_stf);
    return;
  }

  /* Find upper-bound on last output */
  double time_end;
  if (e->policy & engine_policy_cosmology)
    time_end = e->cosmology->a_end * e->delta_time_stf;
  else
    time_end = e->time_end + e->delta_time_stf;

  /* Find next snasphot above current time */
  double time;
  if (e->policy & engine_policy_cosmology)
    time = e->a_first_stf_output;
  else
    time = e->time_first_stf_output;

  int found_stf_time = 0;
  while (time < time_end) {

    /* Output time on the integer timeline */
    if (e->policy & engine_policy_cosmology)
      e->ti_next_stf = log(time / e->cosmology->a_begin) / e->time_base;
    else
      e->ti_next_stf = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (e->ti_next_stf > e->ti_current) {
      found_stf_time = 1;
      break;
    }

    if (e->policy & engine_policy_cosmology)
      time *= e->delta_time_stf;
    else
      time += e->delta_time_stf;
  }

  /* Deal with last snapshot */
  if (!found_stf_time) {
    e->ti_next_stf = -1;
    if (e->verbose) message("No further output time.");
  } else {

    /* Be nice, talk... */
    if (e->policy & engine_policy_cosmology) {
      const float next_stf_time =
          exp(e->ti_next_stf * e->time_base) * e->cosmology->a_begin;
      if (e->verbose)
        message("Next VELOCIraptor time set to a=%e.", next_stf_time);
    } else {
      const float next_stf_time = e->ti_next_stf * e->time_base + e->time_begin;
      if (e->verbose)
        message("Next VELOCIraptor time set to t=%e.", next_stf_time);
    }
  }
}

/**
 * @brief Computes the next time (on the time line) for FoF black holes seeding
 *
 * @param e The #engine.
 */
void engine_compute_next_fof_time(struct engine *e) {

  /* Find upper-bound on last output */
  double time_end;
  if (e->policy & engine_policy_cosmology)
    time_end = e->cosmology->a_end * e->delta_time_fof;
  else
    time_end = e->time_end + e->delta_time_fof;

  /* Find next snasphot above current time */
  double time;
  if (e->policy & engine_policy_cosmology)
    time = e->a_first_fof_call;
  else
    time = e->time_first_fof_call;

  int found_fof_time = 0;
  while (time < time_end) {

    /* Output time on the integer timeline */
    if (e->policy & engine_policy_cosmology)
      e->ti_next_fof = log(time / e->cosmology->a_begin) / e->time_base;
    else
      e->ti_next_fof = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (e->ti_next_fof > e->ti_current) {
      found_fof_time = 1;
      break;
    }

    if (e->policy & engine_policy_cosmology)
      time *= e->delta_time_fof;
    else
      time += e->delta_time_fof;
  }

  /* Deal with last snapshot */
  if (!found_fof_time) {
    e->ti_next_fof = -1;
    if (e->verbose) message("No further FoF time.");
  } else {

    /* Be nice, talk... */
    if (e->policy & engine_policy_cosmology) {
      const float next_fof_time =
          exp(e->ti_next_fof * e->time_base) * e->cosmology->a_begin;
      // if (e->verbose)
      message("Next FoF time set to a=%e.", next_fof_time);
    } else {
      const float next_fof_time = e->ti_next_fof * e->time_base + e->time_begin;
      if (e->verbose) message("Next FoF time set to t=%e.", next_fof_time);
    }
  }
}

/**
 * @brief Initialize all the output_list required by the engine
 *
 * @param e The #engine.
 * @param params The #swift_params.
 */
void engine_init_output_lists(struct engine *e, struct swift_params *params,
                              const struct output_options *output_options) {

  /* Deal with snapshots */
  e->output_list_snapshots = NULL;
  output_list_init(&e->output_list_snapshots, e, "Snapshots",
                   &e->delta_time_snapshot);

  if (e->output_list_snapshots) {

    /* If we are using a different output selection for the
     * various entries, verify that the user did not specify
     * invalid selections. */
    if (e->output_list_snapshots->select_output_on)
      output_list_check_selection(e->output_list_snapshots, output_options);

    engine_compute_next_snapshot_time(e);

    if (e->policy & engine_policy_cosmology)
      e->a_first_snapshot =
          exp(e->ti_next_snapshot * e->time_base) * e->cosmology->a_begin;
    else
      e->time_first_snapshot =
          e->ti_next_snapshot * e->time_base + e->time_begin;
  }

  /* Deal with stats */
  e->output_list_stats = NULL;
  output_list_init(&e->output_list_stats, e, "Statistics",
                   &e->delta_time_statistics);

  if (e->output_list_stats) {
    engine_compute_next_statistics_time(e);

    if (e->policy & engine_policy_cosmology)
      e->a_first_statistics =
          exp(e->ti_next_stats * e->time_base) * e->cosmology->a_begin;
    else
      e->time_first_statistics =
          e->ti_next_stats * e->time_base + e->time_begin;
  }

  /* Deal with stf */
  e->output_list_stf = NULL;
  output_list_init(&e->output_list_stf, e, "StructureFinding",
                   &e->delta_time_stf);

  if (e->output_list_stf) {
    engine_compute_next_stf_time(e);

    if (e->policy & engine_policy_cosmology)
      e->a_first_stf_output =
          exp(e->ti_next_stf * e->time_base) * e->cosmology->a_begin;
    else
      e->time_first_stf_output = e->ti_next_stf * e->time_base + e->time_begin;
  }

  /* Deal with line of sight */
  e->output_list_los = NULL;
  output_list_init(&e->output_list_los, e, "LineOfSight", &e->delta_time_los);

  if (e->output_list_los) {
    engine_compute_next_los_time(e);

    if (e->policy & engine_policy_cosmology)
      e->a_first_los =
          exp(e->ti_next_los * e->time_base) * e->cosmology->a_begin;
    else
      e->time_first_los = e->ti_next_los * e->time_base + e->time_begin;
  }
}
