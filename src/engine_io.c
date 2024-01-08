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

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "active.h"
#include "csds_io.h"
#include "distributed_io.h"
#include "kick.h"
#include "lightcone/lightcone.h"
#include "lightcone/lightcone_array.h"
#include "line_of_sight.h"
#include "parallel_io.h"
#include "power_spectrum.h"
#include "serial_io.h"
#include "single_io.h"
#include "tracers.h"

/* Standard includes */
#include <stdio.h>

/**
 * @brief Finalize the quantities recorded via the snapshot triggers
 *
 * This adds the small amount of time between the particles' last start
 * of step and the current (snapshot) time.
 *
 * @param e The #engine to act on.
 */
void engine_finalize_trigger_recordings(struct engine *e) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);
  struct space *s = e->s;

  /* Finish the recording period for part triggers */
  if (num_snapshot_triggers_part) {
    for (size_t k = 0; k < s->nr_parts; ++k) {

      /* Get a handle on the part. */
      struct part *p = &s->parts[k];
      struct xpart *xp = &s->xparts[k];
      const integertime_t ti_begin =
          get_integer_time_begin(e->ti_current, p->time_bin);

      /* Escape inhibited particles */
      if (part_is_inhibited(p, e)) continue;

      /* We need to escape the special case of a particle that
       * actually ended its time-step on this very step */
      if (e->ti_current - ti_begin == get_integer_timestep(p->time_bin))
        continue;

      /* Time from the start of the particle's step to the snapshot (aka.
       * current time) */
      double missing_time;
      if (with_cosmology) {
        missing_time =
            cosmology_get_delta_time(e->cosmology, ti_begin, e->ti_current);
      } else {
        missing_time = (e->ti_current - ti_begin) * e->time_base;
      }

      tracers_after_timestep_part(
          p, xp, e->internal_units, e->physical_constants, with_cosmology,
          e->cosmology, e->hydro_properties, e->cooling_func, e->time,
          missing_time, e->snapshot_recording_triggers_started_part);
    }
  }

  /* Finish the recording period for spart triggers */
  if (num_snapshot_triggers_spart) {
    for (size_t k = 0; k < s->nr_sparts; ++k) {

      /* Get a handle on the part. */
      struct spart *sp = &s->sparts[k];
      const integertime_t ti_begin =
          get_integer_time_begin(e->ti_current, sp->time_bin);

      /* Escape inhibited particles */
      if (spart_is_inhibited(sp, e)) continue;

      /* We need to escape the special case of a particle that
       * actually ended its time-step on this very step */
      if (e->ti_current - ti_begin == get_integer_timestep(sp->time_bin))
        continue;

      /* Time from the start of the particle's step to the snapshot (aka.
       * current time) */
      double missing_time;
      if (with_cosmology) {
        missing_time =
            cosmology_get_delta_time(e->cosmology, ti_begin, e->ti_current);
      } else {
        missing_time = (e->ti_current - ti_begin) * e->time_base;
      }

      tracers_after_timestep_spart(
          sp, e->internal_units, e->physical_constants, with_cosmology,
          e->cosmology, missing_time,
          e->snapshot_recording_triggers_started_spart);
    }
  }

  /* Finish the recording period for spart triggers */
  if (num_snapshot_triggers_bpart) {
    for (size_t k = 0; k < s->nr_bparts; ++k) {

      /* Get a handle on the part. */
      struct bpart *bp = &s->bparts[k];
      const integertime_t ti_begin =
          get_integer_time_begin(e->ti_current, bp->time_bin);

      /* Escape inhibited particles */
      if (bpart_is_inhibited(bp, e)) continue;

      /* We need to escape the special case of a particle that
       * actually ended its time-step on this very step */
      if (e->ti_current - ti_begin == get_integer_timestep(bp->time_bin))
        continue;

      /* Time from the start of the particle's step to the snapshot (aka.
       * current time) */
      double missing_time;
      if (with_cosmology) {
        missing_time =
            cosmology_get_delta_time(e->cosmology, ti_begin, e->ti_current);
      } else {
        missing_time = (e->ti_current - ti_begin) * e->time_base;
      }

      tracers_after_timestep_bpart(
          bp, e->internal_units, e->physical_constants, with_cosmology,
          e->cosmology, missing_time,
          e->snapshot_recording_triggers_started_bpart);
    }
  }
}

/**
 * @brief dump restart files if it is time to do so and dumps are enabled.
 *
 * @param e the engine.
 * @param drifted_all true if a drift_all has just been performed.
 * @param force force a dump, if dumping is enabled.
 * @return Do we want to stop the run altogether?
 */
int engine_dump_restarts(struct engine *e, const int drifted_all,
                         const int force) {

  /* Are any of the conditions to fully stop a run met? */
  const int end_run_time = e->runtime > e->restart_max_hours_runtime;
  const int stop_file = (e->step % e->restart_stop_steps == 0 &&
                         restart_stop_now(e->restart_dir, 0));

  /* Exit run when told to */
  const int exit_run = (end_run_time || stop_file);

  if (e->restart_dump) {
    ticks tic = getticks();

    const int check_point_time = tic > e->restart_next;

    /* Dump when the time has arrived, or we are told to. */
    int dump = (check_point_time || end_run_time || force || stop_file);

#ifdef WITH_MPI
    /* Synchronize this action from rank 0 (ticks may differ between
     * machines). */
    MPI_Bcast(&dump, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    if (dump) {

      if (e->nodeID == 0) {

        /* Flush the time-step file to avoid gaps in case of crashes
         * before the next automated flush */
        fflush(e->file_timesteps);

        message("Writing restart files");
      }

      /* Clean out the previous saved files, if found. Do this now as we are
       * MPI synchronized. */
      restart_remove_previous(e->restart_file);

      /* Drift all particles first (may have just been done). */
      if (!drifted_all) engine_drift_all(e, /*drift_mpole=*/1);

        /* Free the foreign particles to get more breathing space. */
#ifdef WITH_MPI
      if (e->free_foreign_when_dumping_restart)
        space_free_foreign_parts(e->s, /*clear_cell_pointers=*/1);
#endif

#ifdef WITH_LIGHTCONE
      /* Flush lightcone buffers before dumping restarts */
      lightcone_array_flush(e->lightcone_array_properties, &(e->threadpool),
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

      /* Reallocate freed memory */
      if (e->free_foreign_when_dumping_restart)
        engine_allocate_foreign_particles(e, /*fof=*/0);
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

  /* If we stopped by reaching the time limit, flag that we need to
   * run the resubmission command */
  if (end_run_time && e->resubmit_after_max_hours) e->resubmit = 1;

  return exit_run;
}

/**
 * @brief Writes a snapshot with the current state of the engine
 *
 * @param e The #engine.
 * @param fof Is this a stand-alone FOF call?
 */
void engine_dump_snapshot(struct engine *e, const int fof) {

  return;

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

  /* Finalize the recording periods */
  engine_finalize_trigger_recordings(e);

  /* Get time-step since the last mesh kick */
  if ((e->policy & engine_policy_self_gravity) && e->s->periodic) {
    const int with_cosmology = (e->policy & engine_policy_cosmology);

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

  MPI_Info info;
  MPI_Info_create(&info);

  if (e->snapshot_distributed) {

    write_output_distributed(e, e->internal_units, e->snapshot_units, fof,
                             e->nodeID, e->nr_nodes, MPI_COMM_WORLD, info);

  } else {

#if defined(HAVE_PARALLEL_HDF5)
    write_output_parallel(e, e->internal_units, e->snapshot_units, fof,
                          e->nodeID, e->nr_nodes, MPI_COMM_WORLD, info);
#else
    write_output_serial(e, e->internal_units, e->snapshot_units, fof, e->nodeID,
                        e->nr_nodes, MPI_COMM_WORLD, info);
#endif
  }
  MPI_Info_free(&info);
#else
  write_output_single(e, e->internal_units, e->snapshot_units, fof);
#endif /* WITH_MPI */
#endif /* WITH_HDF5 */

  /* Cancel any triggers that are switched on */
  if (num_snapshot_triggers_part > 0 || num_snapshot_triggers_spart > 0 ||
      num_snapshot_triggers_bpart > 0) {

    /* Reset the trigger flags */
    for (int i = 0; i < num_snapshot_triggers_part; ++i)
      e->snapshot_recording_triggers_started_part[i] = 0;
    for (int i = 0; i < num_snapshot_triggers_spart; ++i)
      e->snapshot_recording_triggers_started_spart[i] = 0;
    for (int i = 0; i < num_snapshot_triggers_bpart; ++i)
      e->snapshot_recording_triggers_started_bpart[i] = 0;

    /* Reser the tracers themselves */
    space_after_snap_tracer(e->s, e->verbose);
  }

  /* Flag that we dumped a snapshot */
  e->step_props |= engine_step_prop_snapshot;

  clocks_gettime(&time2);
  if (e->verbose)
    message("writing particle properties took %.3f %s.",
            (float)clocks_diff(&time1, &time2), clocks_getunit());

  /* Run the post-dump command if required */
  if (e->nodeID == 0) {
    engine_run_on_dump(e);
  }
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
void engine_io(struct engine *e) {
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int with_stf = (e->policy & engine_policy_structure_finding);
  const int with_los = (e->policy & engine_policy_line_of_sight);
  const int with_fof = (e->policy & engine_policy_fof);
  const int with_power = (e->policy & engine_policy_power_spectra);

  /* What kind of output are we getting? */
  enum output_type {
    output_none,
    output_snapshot,
    output_statistics,
    output_ps,
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

  /* Do we want a power-spectrum? */
  if (e->ti_end_min > e->ti_next_ps && e->ti_next_ps > 0) {
    if (e->ti_next_ps < ti_output) {
      ti_output = e->ti_next_ps;
      type = output_ps;
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

          /* Free the foreign particles to get more breathing space.
           * If called, the FOF code itself will reallocate what it needs. */
#ifdef WITH_MPI
        space_free_foreign_parts(e->s, /*clear_cell_pointers=*/1);
#endif

        /* Do we want FoF group IDs in the snapshot? */
        if (with_fof && e->snapshot_invoke_fof) {
          engine_fof(e, /*dump_results=*/1, /*dump_debug=*/0,
                     /*seed_black_holes=*/0, /*buffers allocated=*/0);
        }

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

        /* Do we want power spectrum outputs? */
        if (with_power && e->snapshot_invoke_ps) {
          calc_all_power_spectra(e->power_data, e->s, &e->threadpool,
                                 e->verbose);
        }

        /* Dump... */
        engine_dump_snapshot(e, /*fof=*/0);

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
        engine_allocate_foreign_particles(e, /*fof=*/0);
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
        engine_allocate_foreign_particles(e, /*fof=*/0);
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

      case output_ps:

        /* Compute the PS */
        calc_all_power_spectra(e->power_data, e->s, &e->threadpool, e->verbose);

        /* Move on */
        engine_compute_next_ps_time(e);

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

    /* Do we want a power-spectrum? */
    if (e->ti_end_min > e->ti_next_ps && e->ti_next_ps > 0) {
      if (e->ti_next_ps < ti_output) {
        ti_output = e->ti_next_ps;
        type = output_ps;
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

    /* Time until the next snapshot */
    double time_to_next_snap;
    if (e->policy & engine_policy_cosmology) {
      time_to_next_snap = cosmology_get_delta_time(e->cosmology, e->ti_current,
                                                   e->ti_next_snapshot);
    } else {
      time_to_next_snap = (e->ti_next_snapshot - e->ti_current) * e->time_base;
    }

    /* Do we need to reduce any of the recording trigger times? */
    for (int k = 0; k < num_snapshot_triggers_part; ++k) {
      if (e->snapshot_recording_triggers_desired_part[k] > 0) {
        if (e->snapshot_recording_triggers_desired_part[k] >
            time_to_next_snap) {
          e->snapshot_recording_triggers_part[k] = time_to_next_snap;
        } else {
          e->snapshot_recording_triggers_part[k] =
              e->snapshot_recording_triggers_desired_part[k];
        }
      }
    }
    for (int k = 0; k < num_snapshot_triggers_spart; ++k) {
      if (e->snapshot_recording_triggers_desired_spart[k] > 0) {
        if (e->snapshot_recording_triggers_desired_spart[k] >
            time_to_next_snap) {
          e->snapshot_recording_triggers_spart[k] = time_to_next_snap;
        } else {
          e->snapshot_recording_triggers_spart[k] =
              e->snapshot_recording_triggers_desired_spart[k];
        }
      }
    }
    for (int k = 0; k < num_snapshot_triggers_bpart; ++k) {
      if (e->snapshot_recording_triggers_desired_bpart[k] > 0) {
        if (e->snapshot_recording_triggers_desired_bpart[k] >
            time_to_next_snap) {
          e->snapshot_recording_triggers_bpart[k] = time_to_next_snap;
        } else {
          e->snapshot_recording_triggers_bpart[k] =
              e->snapshot_recording_triggers_desired_bpart[k];
        }
      }
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
      if (e->verbose) message("Next FoF time set to a=%e.", next_fof_time);
    } else {
      const float next_fof_time = e->ti_next_fof * e->time_base + e->time_begin;
      if (e->verbose) message("Next FoF time set to t=%e.", next_fof_time);
    }
  }
}

/**
 * @brief Computes the next time (on the time line) for a power-spectrum dump
 *
 * @param e The #engine.
 */
void engine_compute_next_ps_time(struct engine *e) {
  /* Do output_list file case */
  if (e->output_list_ps) {
    output_list_read_next_time(e->output_list_ps, e, "power spectrum",
                               &e->ti_next_ps);
    return;
  }

  /* Find upper-bound on last output */
  double time_end;
  if (e->policy & engine_policy_cosmology)
    time_end = e->cosmology->a_end * e->delta_time_ps;
  else
    time_end = e->time_end + e->delta_time_ps;

  /* Find next ps above current time */
  double time;
  if (e->policy & engine_policy_cosmology)
    time = e->a_first_ps_output;
  else
    time = e->time_first_ps_output;

  int found_ps_time = 0;
  while (time < time_end) {

    /* Output time on the integer timeline */
    if (e->policy & engine_policy_cosmology)
      e->ti_next_ps = log(time / e->cosmology->a_begin) / e->time_base;
    else
      e->ti_next_ps = (time - e->time_begin) / e->time_base;

    /* Found it? */
    if (e->ti_next_ps > e->ti_current) {
      found_ps_time = 1;
      break;
    }

    if (e->policy & engine_policy_cosmology)
      time *= e->delta_time_ps;
    else
      time += e->delta_time_ps;
  }

  /* Deal with last line of sight */
  if (!found_ps_time) {
    e->ti_next_ps = -1;
    if (e->verbose) message("No further PS output time.");
  } else {

    /* Be nice, talk... */
    if (e->policy & engine_policy_cosmology) {
      const double next_ps_time =
          exp(e->ti_next_ps * e->time_base) * e->cosmology->a_begin;
      if (e->verbose)
        message("Next output time for power spectrum set to a=%e.",
                next_ps_time);
    } else {
      const double next_ps_time = e->ti_next_ps * e->time_base + e->time_begin;
      if (e->verbose)
        message("Next output time for power spectrum set to t=%e.",
                next_ps_time);
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
  if (e->policy & engine_policy_structure_finding) {

    e->output_list_stf = NULL;
    output_list_init(&e->output_list_stf, e, "StructureFinding",
                     &e->delta_time_stf);

    if (e->output_list_stf) {
      engine_compute_next_stf_time(e);

      if (e->policy & engine_policy_cosmology)
        e->a_first_stf_output =
            exp(e->ti_next_stf * e->time_base) * e->cosmology->a_begin;
      else
        e->time_first_stf_output =
            e->ti_next_stf * e->time_base + e->time_begin;
    }
  }

  /* Deal with line of sight */
  if (e->policy & engine_policy_line_of_sight) {

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

  /* Deal with power-spectra */
  if (e->policy & engine_policy_power_spectra) {

    e->output_list_ps = NULL;
    output_list_init(&e->output_list_ps, e, "PowerSpectrum", &e->delta_time_ps);

    if (e->output_list_ps) {
      engine_compute_next_ps_time(e);

      if (e->policy & engine_policy_cosmology)
        e->a_first_ps_output =
            exp(e->ti_next_ps * e->time_base) * e->cosmology->a_begin;
      else
        e->time_first_ps_output = e->ti_next_ps * e->time_base + e->time_begin;
    }
  }
}

/**
 * @brief Checks whether we passed a certain delta time before the next
 * snapshot and need to trigger a recording.
 *
 * If a recording has to start, we also loop over the particles to correct
 * for the time between the particles' end of time-step and the actual start
 * of trigger.
 *
 * @param e The #engine.
 */
void engine_io_check_snapshot_triggers(struct engine *e) {

  struct space *s = e->s;
  const int with_cosmology = (e->policy & engine_policy_cosmology);

  /* Time until the next snapshot */
  double time_to_next_snap;
  if (e->policy & engine_policy_cosmology) {
    time_to_next_snap = cosmology_get_delta_time(e->cosmology, e->ti_current,
                                                 e->ti_next_snapshot);
  } else {
    time_to_next_snap = (e->ti_next_snapshot - e->ti_current) * e->time_base;
  }

  /* Should any not yet switched on trigger be activated? (part version) */
  for (int i = 0; i < num_snapshot_triggers_part; ++i) {

    if (time_to_next_snap <= e->snapshot_recording_triggers_part[i] &&
        e->snapshot_recording_triggers_part[i] > 0. &&
        !e->snapshot_recording_triggers_started_part[i]) {
      e->snapshot_recording_triggers_started_part[i] = 1;

      /* Be vocal about this */
      if (e->verbose)
        message(
            "Snapshot will be dumped in %e U_t. Recording trigger for part "
            "activated.",
            e->snapshot_recording_triggers_part[i]);

      /* We now need to loop over the particles to preemptively deduct the
       * extra time logged between the particles' start of step and the
       * actual start of the trigger */
      for (size_t k = 0; k < s->nr_parts; ++k) {

        /* Get a handle on the part. */
        struct part *p = &s->parts[k];
        struct xpart *xp = &s->xparts[k];
        const integertime_t ti_begin =
            get_integer_time_begin(e->ti_current, p->time_bin);

        /* Escape inhibited particles */
        if (part_is_inhibited(p, e)) continue;

        /* Time from the start of the particle's step to the snapshot */
        double total_time;
        if (with_cosmology) {
          total_time = cosmology_get_delta_time(e->cosmology, ti_begin,
                                                e->ti_next_snapshot);
        } else {
          total_time = (e->ti_next_snapshot - ti_begin) * e->time_base;
        }

        /* Time to deduct = time since the start of the step - trigger time */
        const double time_to_remove =
            total_time - e->snapshot_recording_triggers_part[i];

#ifdef SWIFT_DEBUG_CHECKS
        if (time_to_remove < 0.)
          error("Invalid time to deduct! %e", time_to_remove);
#endif

        /* Note that we need to use a separate array (not the raw
         * e->snapshot_recording_triggers_part) as we only want to
         * update one entry */
        int my_temp_array[num_snapshot_triggers_part];
        memset(my_temp_array, 0, sizeof(int) * num_snapshot_triggers_part);
        my_temp_array[i] = 1;

        tracers_after_timestep_part(
            p, xp, e->internal_units, e->physical_constants, with_cosmology,
            e->cosmology, e->hydro_properties, e->cooling_func, e->time,
            -time_to_remove, my_temp_array);
      }
    }
  }

  /* Should any not yet switched on trigger be activated? (spart version) */
  for (int i = 0; i < num_snapshot_triggers_spart; ++i) {

    if (time_to_next_snap <= e->snapshot_recording_triggers_spart[i] &&
        e->snapshot_recording_triggers_spart[i] > 0. &&
        !e->snapshot_recording_triggers_started_spart[i]) {
      e->snapshot_recording_triggers_started_spart[i] = 1;

      /* Be vocal about this */
      if (e->verbose)
        message(
            "Snapshot will be dumped in %e U_t. Recording trigger for spart "
            "activated.",
            e->snapshot_recording_triggers_spart[i]);

      /* We now need to loop over the particles to preemptively deduct the
       * extra time logged between the particles' start of step and the
       * actual start of the trigger */
      for (size_t k = 0; k < s->nr_sparts; ++k) {

        /* Get a handle on the spart. */
        struct spart *sp = &s->sparts[k];
        const integertime_t ti_begin =
            get_integer_time_begin(e->ti_current, sp->time_bin);

        /* Escape inhibited particles */
        if (spart_is_inhibited(sp, e)) continue;

        /* Time from the start of the particle's step to the snapshot */
        double total_time;
        if (with_cosmology) {
          total_time = cosmology_get_delta_time(e->cosmology, ti_begin,
                                                e->ti_next_snapshot);
        } else {
          total_time = (e->ti_next_snapshot - ti_begin) * e->time_base;
        }

        /* Time to deduct = time since the start of the step - trigger time */
        const double time_to_remove =
            total_time - e->snapshot_recording_triggers_part[i];

#ifdef SWIFT_DEBUG_CHECKS
        if (time_to_remove < 0.)
          error("Invalid time to deduct! %e", time_to_remove);
#endif

        /* Note that we need to use a separate array (not the raw
         * e->snapshot_recording_triggers_part) as we only want to
         * update one entry */
        int my_temp_array[num_snapshot_triggers_spart];
        memset(my_temp_array, 0, sizeof(int) * num_snapshot_triggers_spart);
        my_temp_array[i] = 1;

        tracers_after_timestep_spart(
            sp, e->internal_units, e->physical_constants, with_cosmology,
            e->cosmology, -time_to_remove, my_temp_array);
      }
    }
  }

  /* Should any not yet switched on trigger be activated? (bpart version) */
  for (int i = 0; i < num_snapshot_triggers_bpart; ++i) {

    if (time_to_next_snap <= e->snapshot_recording_triggers_bpart[i] &&
        e->snapshot_recording_triggers_bpart[i] > 0. &&
        !e->snapshot_recording_triggers_started_bpart[i]) {
      e->snapshot_recording_triggers_started_bpart[i] = 1;

      /* Be vocal about this */
      if (e->verbose)
        message(
            "Snapshot will be dumped in %e U_t. Recording trigger for bpart "
            "activated.",
            e->snapshot_recording_triggers_bpart[i]);

      /* We now need to loop over the particles to preemptively deduct the
       * extra time logged between the particles' start of step and the
       * actual start of the trigger */
      for (size_t k = 0; k < s->nr_bparts; ++k) {

        /* Get a handle on the bpart. */
        struct bpart *bp = &s->bparts[k];
        const integertime_t ti_begin =
            get_integer_time_begin(e->ti_current, bp->time_bin);

        /* Escape inhibited particles */
        if (bpart_is_inhibited(bp, e)) continue;

        /* Time from the start of the particle's step to the snapshot */
        double total_time;
        if (with_cosmology) {
          total_time = cosmology_get_delta_time(e->cosmology, ti_begin,
                                                e->ti_next_snapshot);
        } else {
          total_time = (e->ti_next_snapshot - ti_begin) * e->time_base;
        }

        /* Time to deduct = time since the start of the step - trigger time */
        const double time_to_remove =
            total_time - e->snapshot_recording_triggers_part[i];

#ifdef SWIFT_DEBUG_CHECKS
        if (time_to_remove < 0.)
          error("Invalid time to deduct! %e", time_to_remove);
#endif

        /* Note that we need to use a separate array (not the raw
         * e->snapshot_recording_triggers_part) as we only want to
         * update one entry */
        int my_temp_array[num_snapshot_triggers_bpart];
        memset(my_temp_array, 0, sizeof(int) * num_snapshot_triggers_bpart);
        my_temp_array[i] = 1;

        tracers_after_timestep_bpart(
            bp, e->internal_units, e->physical_constants, with_cosmology,
            e->cosmology, -time_to_remove, my_temp_array);
      }
    }
  }
}
