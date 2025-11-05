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

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "active.h"
#include "lightcone/lightcone_array.h"
#include "star_formation_logger.h"
#include "timeline.h"

/**
 * @brief Data collected from the cells at the end of a time-step
 */
struct end_of_step_data {

  size_t updated, g_updated, s_updated, sink_updated, b_updated;
  size_t inhibited, g_inhibited, s_inhibited, sink_inhibited, b_inhibited;
  integertime_t ti_hydro_end_min, ti_hydro_beg_max;
  integertime_t ti_rt_end_min, ti_rt_beg_max;
  integertime_t ti_gravity_end_min, ti_gravity_beg_max;
  integertime_t ti_stars_end_min, ti_stars_beg_max;
  integertime_t ti_sinks_end_min, ti_sinks_beg_max;
  integertime_t ti_black_holes_end_min, ti_black_holes_beg_max;
  struct engine *e;
  struct star_formation_history sfh;
  float runtime;
  int flush_lightcone_maps;
  double deadtime;
  float csds_file_size_gb;
};

/**
 * @brief Mapping function to collect the data from the end of the step
 *
 * This function will call a recursive function on all the top-level cells
 * to collect the information we are after.
 *
 * @param map_data The list of cells with tasks on this node.
 * @param num_elements The number of elements in the list this thread will work
 * on.
 * @param extra_data The #engine.
 */
void engine_collect_end_of_step_mapper(void *map_data, int num_elements,
                                       void *extra_data) {

  struct end_of_step_data *data = (struct end_of_step_data *)extra_data;
  const struct engine *e = data->e;
  struct space *s = e->s;
  int *local_cells = (int *)map_data;
  struct star_formation_history *sfh_top = &data->sfh;

  /* Local collectible */
  size_t updated = 0, g_updated = 0, s_updated = 0, sink_updated = 0,
         b_updated = 0;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_beg_max = 0;
  integertime_t ti_rt_end_min = max_nr_timesteps, ti_rt_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps, ti_stars_beg_max = 0;
  integertime_t ti_sinks_end_min = max_nr_timesteps, ti_sinks_beg_max = 0;
  integertime_t ti_black_holes_end_min = max_nr_timesteps,
                ti_black_holes_beg_max = 0;

  /* Local Star formation history properties */
  struct star_formation_history sfh_updated;

  /* Initialize the star formation structs for this engine to zero */
  star_formation_logger_init(&sfh_updated);

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &s->cells_top[local_cells[ind]];

    if (c->hydro.count > 0 || c->grav.count > 0 || c->stars.count > 0 ||
        c->black_holes.count > 0 || c->sinks.count > 0) {

      /* Aggregate data */
      if (c->hydro.ti_end_min > e->ti_current)
        ti_hydro_end_min = min(ti_hydro_end_min, c->hydro.ti_end_min);
      ti_hydro_beg_max = max(ti_hydro_beg_max, c->hydro.ti_beg_max);

      if (c->rt.ti_rt_end_min > e->ti_current)
        ti_rt_end_min = min(c->rt.ti_rt_end_min, ti_rt_end_min);
      ti_rt_beg_max = max(c->rt.ti_rt_beg_max, ti_rt_beg_max);

      if (c->grav.ti_end_min > e->ti_current)
        ti_gravity_end_min = min(ti_gravity_end_min, c->grav.ti_end_min);
      ti_gravity_beg_max = max(ti_gravity_beg_max, c->grav.ti_beg_max);

      if (c->stars.ti_end_min > e->ti_current)
        ti_stars_end_min = min(ti_stars_end_min, c->stars.ti_end_min);
      ti_stars_beg_max = max(ti_stars_beg_max, c->stars.ti_beg_max);

      if (c->sinks.ti_end_min > e->ti_current)
        ti_sinks_end_min = min(ti_sinks_end_min, c->sinks.ti_end_min);
      ti_sinks_beg_max = max(ti_sinks_beg_max, c->sinks.ti_beg_max);

      if (c->black_holes.ti_end_min > e->ti_current)
        ti_black_holes_end_min =
            min(ti_black_holes_end_min, c->black_holes.ti_end_min);
      ti_black_holes_beg_max =
          max(ti_black_holes_beg_max, c->black_holes.ti_beg_max);

      updated += c->hydro.updated;
      g_updated += c->grav.updated;
      s_updated += c->stars.updated;
      sink_updated += c->sinks.updated;
      b_updated += c->black_holes.updated;

      /* Check if the cell was inactive and in that case reorder the SFH */
      if (!cell_is_starting_hydro(c, e)) {
        star_formation_logger_log_inactive_cell(&c->stars.sfh);
      }

      /* Get the star formation history from the current cell and store it in
       * the star formation history struct */
      star_formation_logger_add(&sfh_updated, &c->stars.sfh);

      /* Collected, so clear for next time. */
      c->hydro.updated = 0;
      c->grav.updated = 0;
      c->stars.updated = 0;
      c->sinks.updated = 0;
      c->black_holes.updated = 0;
    }
  }

  /* Let's write back to the global data.
   * We use the space lock to garanty single access*/
  if (lock_lock(&s->lock) == 0) {
    data->updated += updated;
    data->g_updated += g_updated;
    data->s_updated += s_updated;
    data->sink_updated += sink_updated;
    data->b_updated += b_updated;

    /* Add the SFH information from this engine to the global data */
    star_formation_logger_add(sfh_top, &sfh_updated);

    if (ti_hydro_end_min > e->ti_current)
      data->ti_hydro_end_min = min(ti_hydro_end_min, data->ti_hydro_end_min);
    data->ti_hydro_beg_max = max(ti_hydro_beg_max, data->ti_hydro_beg_max);

    if (ti_rt_end_min > e->ti_current)
      data->ti_rt_end_min = min(ti_rt_end_min, data->ti_rt_end_min);
    data->ti_rt_beg_max = max(ti_rt_beg_max, data->ti_rt_beg_max);

    if (ti_gravity_end_min > e->ti_current)
      data->ti_gravity_end_min =
          min(ti_gravity_end_min, data->ti_gravity_end_min);
    data->ti_gravity_beg_max =
        max(ti_gravity_beg_max, data->ti_gravity_beg_max);

    if (ti_stars_end_min > e->ti_current)
      data->ti_stars_end_min = min(ti_stars_end_min, data->ti_stars_end_min);
    data->ti_stars_beg_max = max(ti_stars_beg_max, data->ti_stars_beg_max);

    if (ti_sinks_end_min > e->ti_current)
      data->ti_sinks_end_min = min(ti_sinks_end_min, data->ti_sinks_end_min);
    data->ti_sinks_beg_max = max(ti_sinks_beg_max, data->ti_sinks_beg_max);

    if (ti_black_holes_end_min > e->ti_current)
      data->ti_black_holes_end_min =
          min(ti_black_holes_end_min, data->ti_black_holes_end_min);
    data->ti_black_holes_beg_max =
        max(ti_black_holes_beg_max, data->ti_black_holes_beg_max);
  }

  if (lock_unlock(&s->lock) != 0) error("Failed to unlock the space");
}

/**
 * @brief Collects the next time-step and rebuild flag.
 *
 * The next time-step is determined by making each super-cell recurse to
 * collect the minimal of ti_end and the number of updated particles.  When in
 * MPI mode this routines reduces these across all nodes and also collects the
 * forcerebuild flag -- this is so that we only use a single collective MPI
 * call per step for all these values.
 *
 * Note that the results are stored in e->collect_group1 struct not in the
 * engine fields, unless apply is true. These can be applied field-by-field
 * or all at once using collectgroup1_copy();
 *
 * @param e The #engine.
 * @param apply whether to apply the results to the engine or just keep in the
 *              group1 struct.
 */
void engine_collect_end_of_step(struct engine *e, int apply) {

  const ticks tic = getticks();
  struct space *s = e->s;
  struct end_of_step_data data;
  data.updated = 0, data.g_updated = 0, data.s_updated = 0, data.b_updated = 0;
  data.sink_updated = 0;
  data.ti_hydro_end_min = max_nr_timesteps, data.ti_hydro_beg_max = 0;
  data.ti_rt_end_min = max_nr_timesteps, data.ti_rt_beg_max = 0;
  data.ti_gravity_end_min = max_nr_timesteps, data.ti_gravity_beg_max = 0;
  data.ti_stars_end_min = max_nr_timesteps, data.ti_stars_beg_max = 0;
  data.ti_sinks_end_min = max_nr_timesteps, data.ti_sinks_beg_max = 0;
  data.ti_black_holes_end_min = max_nr_timesteps,
  data.ti_black_holes_beg_max = 0, data.e = e, data.csds_file_size_gb = 0;

#ifdef WITH_CSDS
  /* Get the file size from the CSDS. */
  if (e->policy & engine_policy_csds)
    data.csds_file_size_gb =
        csds_logfile_writer_get_current_filesize_used_gb(&e->csds->logfile);
#endif

  /* Need to use a consistent check of the hours since we started. */
  data.runtime = clocks_get_hours_since_start();

  /* Get flag to determine if lightcone maps buffers should be flushed on this
   * step */
  data.flush_lightcone_maps =
      lightcone_array_trigger_map_update(e->lightcone_array_properties);

  data.deadtime = e->local_deadtime;

  /* Initialize the total SFH of the simulation to zero */
  star_formation_logger_init(&data.sfh);

  /* Collect information from the local top-level cells */
  threadpool_map(&e->threadpool, engine_collect_end_of_step_mapper,
                 s->local_cells_top, s->nr_local_cells, sizeof(int),
                 threadpool_auto_chunk_size, &data);

  /* Get the number of inhibited particles from the space-wide counters
   * since these have been updated atomically during the time-steps. */
  data.inhibited = s->nr_inhibited_parts;
  data.g_inhibited = s->nr_inhibited_gparts;
  data.s_inhibited = s->nr_inhibited_sparts;
  data.sink_inhibited = s->nr_inhibited_sinks;
  data.b_inhibited = s->nr_inhibited_bparts;

  /* Store these in the temporary collection group. */
  collectgroup1_init(
      &e->collect_group1, data.updated, data.g_updated, data.s_updated,
      data.sink_updated, data.b_updated, data.inhibited, data.g_inhibited,
      data.s_inhibited, data.sink_inhibited, data.b_inhibited,
      data.ti_hydro_end_min, data.ti_hydro_beg_max, data.ti_rt_end_min,
      data.ti_rt_beg_max, data.ti_gravity_end_min, data.ti_gravity_beg_max,
      data.ti_stars_end_min, data.ti_stars_beg_max, data.ti_sinks_end_min,
      data.ti_sinks_beg_max, data.ti_black_holes_end_min,
      data.ti_black_holes_beg_max, e->forcerebuild, e->s->tot_cells,
      e->sched.nr_tasks, (float)e->sched.nr_tasks / (float)e->s->tot_cells,
      data.sfh, data.runtime, data.flush_lightcone_maps, data.deadtime,
      data.csds_file_size_gb);

/* Aggregate collective data from the different nodes for this step. */
#ifdef WITH_MPI
  collectgroup1_reduce(&e->collect_group1);

#ifdef SWIFT_DEBUG_CHECKS
  {
    /* Check the above using the original MPI calls. */
    integertime_t in_i[2], out_i[2];
    in_i[0] = 0;
    in_i[1] = 0;
    out_i[0] = data.ti_hydro_end_min;
    out_i[1] = data.ti_gravity_end_min;
    if (MPI_Allreduce(out_i, in_i, 2, MPI_LONG_LONG_INT, MPI_MIN,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate ti_end_min.");
    if (in_i[0] != (long long)e->collect_group1.ti_hydro_end_min)
      error("Failed to get same ti_hydro_end_min, is %lld, should be %lld",
            in_i[0], e->collect_group1.ti_hydro_end_min);
    if (in_i[1] != (long long)e->collect_group1.ti_gravity_end_min)
      error("Failed to get same ti_gravity_end_min, is %lld, should be %lld",
            in_i[1], e->collect_group1.ti_gravity_end_min);

    long long in_ll[4], out_ll[4];
    out_ll[0] = data.updated;
    out_ll[1] = data.g_updated;
    out_ll[2] = data.s_updated;
    out_ll[3] = data.b_updated;
    if (MPI_Allreduce(out_ll, in_ll, 4, MPI_LONG_LONG_INT, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate particle counts.");
    if (in_ll[0] != (long long)e->collect_group1.updated)
      error("Failed to get same updated, is %lld, should be %lld", in_ll[0],
            e->collect_group1.updated);
    if (in_ll[1] != (long long)e->collect_group1.g_updated)
      error("Failed to get same g_updated, is %lld, should be %lld", in_ll[1],
            e->collect_group1.g_updated);
    if (in_ll[2] != (long long)e->collect_group1.s_updated)
      error("Failed to get same s_updated, is %lld, should be %lld", in_ll[2],
            e->collect_group1.s_updated);
    if (in_ll[3] != (long long)e->collect_group1.b_updated)
      error("Failed to get same b_updated, is %lld, should be %lld", in_ll[3],
            e->collect_group1.b_updated);

    out_ll[0] = data.inhibited;
    out_ll[1] = data.g_inhibited;
    out_ll[2] = data.s_inhibited;
    out_ll[3] = data.b_inhibited;
    if (MPI_Allreduce(out_ll, in_ll, 4, MPI_LONG_LONG_INT, MPI_SUM,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate particle counts.");
    if (in_ll[0] != (long long)e->collect_group1.inhibited)
      error("Failed to get same inhibited, is %lld, should be %lld", in_ll[0],
            e->collect_group1.inhibited);
    if (in_ll[1] != (long long)e->collect_group1.g_inhibited)
      error("Failed to get same g_inhibited, is %lld, should be %lld", in_ll[1],
            e->collect_group1.g_inhibited);
    if (in_ll[2] != (long long)e->collect_group1.s_inhibited)
      error("Failed to get same s_inhibited, is %lld, should be %lld", in_ll[2],
            e->collect_group1.s_inhibited);
    if (in_ll[3] != (long long)e->collect_group1.b_inhibited)
      error("Failed to get same b_inhibited, is %lld, should be %lld", in_ll[3],
            e->collect_group1.b_inhibited);

    int buff = 0;
    if (MPI_Allreduce(&e->forcerebuild, &buff, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD) != MPI_SUCCESS)
      error("Failed to aggregate the rebuild flag across nodes.");
    if (!!buff != !!e->collect_group1.forcerebuild)
      error(
          "Failed to get same rebuild flag from all nodes, is %d,"
          "should be %d",
          buff, e->collect_group1.forcerebuild);
  }
#endif
#endif

  /* Apply to the engine, if requested. */
  if (apply) collectgroup1_apply(&e->collect_group1, e);

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Mapping function to collect the data from the end of the sub-cycle
 *
 * @param map_data The list of cells with tasks on this node.
 * @param num_elements The number of elements in the list this thread will work
 * on.
 * @param extra_data The #engine.
 */
void engine_collect_end_of_sub_cycle_mapper(void *map_data, int num_elements,
                                            void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  int *local_cells = (int *)map_data;

  /* Local collectible */
  long long rt_updated = 0LL;

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &s->cells_top[local_cells[ind]];

    if (c->hydro.count > 0) {

      /* Aggregate data */
      rt_updated += c->rt.updated;

      /* Collected, so clear for next time. */
      c->rt.updated = 0;
    }
  }

  /* write back to the global data. */
  atomic_add(&e->rt_updates, rt_updated);
}

/**
 * @brief Collects additional data at the end of a subcycle.
 * This function does not collect any data relevant to the
 * time-steps or time integration.
 *
 * @param e The #engine.
 */
void engine_collect_end_of_sub_cycle(struct engine *e) {

  const ticks tic = getticks();
  struct space *s = e->s;

  /* Collect information from the local top-level cells */
  threadpool_map(&e->threadpool, engine_collect_end_of_sub_cycle_mapper,
                 s->local_cells_top, s->nr_local_cells, sizeof(int),
                 threadpool_auto_chunk_size, e);

#ifdef WITH_MPI

  /* Aggregate collective data from the different nodes for this step. */
  int test;
  long long rt_updates_tot = 0ll;
  test = MPI_Reduce(&e->rt_updates, &rt_updates_tot, 1, MPI_LONG_LONG, MPI_SUM,
                    0, MPI_COMM_WORLD);
  if (test != MPI_SUCCESS) error("MPI reduce failed");

  double global_deadtime = 0.;
  test = MPI_Reduce(&e->local_deadtime, &global_deadtime, 1, MPI_DOUBLE,
                    MPI_SUM, 0, MPI_COMM_WORLD);
  if (test != MPI_SUCCESS) error("MPI reduce failed");

  /* Overwrite only on rank 0. */
  if (e->nodeID == 0) {
    e->rt_updates = rt_updates_tot;
    e->global_deadtime = global_deadtime;
  }

#else

  e->global_deadtime = e->local_deadtime;

#endif

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
